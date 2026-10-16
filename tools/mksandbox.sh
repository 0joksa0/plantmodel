#!/usr/bin/env bash
# Create an isolated copy of the repository for one experiment.
#
# Why this exists: attribution requires that each change be measured on its own,
# from a known starting point. Editing src/ in place and running a sweep gives a
# number that cannot be attributed to anything, because the previous experiment's
# edits are still there. This produces a clean tree at a chosen commit, which is
# deleted when the experiment is done.
#
# `git archive` rather than `git worktree`: a worktree registers metadata under
# .git/worktrees, and in a sandboxed environment where the repository is a bind
# mount those files cannot be removed afterwards ("Device or resource busy"),
# leaving stale registrations that block the name. An archive export is a plain
# directory -- no git metadata, removable with rm -rf, and identical in content.
# It costs ~22 MB per copy.
#
# The trap this script exists to close: *.csv is in .gitignore, so an export has
# an EMPTY data/paper/ and the runners then read defaults instead of the
# reference series -- a wrong model with no error message. The gitignored inputs
# are copied in explicitly below.
#
# Usage:
#   tools/mksandbox.sh <name> [commit-ish]
#   tools/mksandbox.sh baseline HEAD
#   tools/mksandbox.sh kappa_only 073f4ed
#
# Removal:
#   rm -rf .sandboxes/<name>
#
set -euo pipefail

NAME="${1:?usage: mksandbox.sh <name> [commit-ish]}"
REF="${2:-HEAD}"

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"

DEST=".sandboxes/${NAME}"

if [ -e "$DEST" ]; then
    echo "error: $DEST already exists. Remove it first:" >&2
    echo "  rm -rf $DEST" >&2
    exit 1
fi

mkdir -p "$DEST"
git archive --format=tar "$REF" | tar -x -C "$DEST"

# Gitignored inputs the runners read at startup. Without these an export looks
# complete and produces wrong numbers rather than an error.
#
# MEASUREMENT DATA is copied from the working tree -- it is not tracked, so the
# working tree is the only copy that exists.
SEED=(
    "data/paper/sugars_starch_reference.csv"
    "data/paper/N_data.csv"
    "data/paper/nitrogen_response_reference.csv"
)
for f in "${SEED[@]}"; do
    if [ -f "$f" ]; then
        mkdir -p "$DEST/$(dirname "$f")"
        cp "$f" "$DEST/$f"
    else
        echo "warning: seed file missing in main tree: $f" >&2
    fi
done

# CONFIGURATION comes from the COMMIT, never from the working tree.
#
# This distinction was got wrong once and it invalidated an anchor measurement.
# config/default.config is gitignored, so an earlier version of this script
# copied it from the working tree "because that is the only copy" -- but the
# working tree carried an uncommitted parameter edit
# (sucrose_consumption_transport 0.0035 -> 0.035). Every sandbox therefore ran
# that edit while claiming to be clean HEAD, and the effect being attributed was
# already present in the reference. `git show <ref>:<path>` reads the file as the
# commit has it, which is what "isolated copy at commit X" has to mean.
for f in config/default.config config/sensitivity_params.csv config/sweep_params.csv; do
    if git cat-file -e "${REF}:${f}" 2>/dev/null; then
        mkdir -p "$DEST/$(dirname "$f")"
        git show "${REF}:${f}" > "$DEST/$f"
    elif [ -f "$f" ]; then
        # Not in the commit at all: fall back to the working tree, and SAY SO,
        # because the sandbox is then not a pure export.
        mkdir -p "$DEST/$(dirname "$f")"
        cp "$f" "$DEST/$f"
        echo "note: $f is untracked at $REF; seeded from the working tree" >&2
    fi
done

# Per-photoperiod ablation configs, also gitignored in part.
if [ -d config/ablation ]; then
    mkdir -p "$DEST/config"
    cp -r config/ablation "$DEST/config/"
fi

# The yardstick. Every sandbox is measured with the same tool, so it is copied in
# and registered against the worktree's OWN CMakeLists -- not the main tree's,
# which references untracked targets that do not exist at this commit.
if [ -f tools/model_report.c ]; then
    cp tools/model_report.c "$DEST/tools/"
    if ! grep -q "model_report" "$DEST/CMakeLists.txt"; then
        # Feature probes. The yardstick must build against an UNMODIFIED tree as
        # well as against variants that add fields, so optional refinements are
        # guarded by these macros rather than assumed present.
        FEAT=""
        grep -q "leaf_temperature_day"  "$DEST/include/model/input.h" 2>/dev/null \
            && FEAT="$FEAT -DPLANTMODEL_HAS_DAYNIGHT_TEMP"
        grep -q "light_profile"         "$DEST/include/model/input.h" 2>/dev/null \
            && FEAT="$FEAT -DPLANTMODEL_HAS_LIGHT_PROFILE"
        grep -q "diverged"              "$DEST/include/simulation/simulation.h" 2>/dev/null \
            && FEAT="$FEAT -DPLANTMODEL_HAS_DIVERGED_FLAG"

        cat >> "$DEST/CMakeLists.txt" <<CMAKE

add_executable(model_report \${SRC} tools/model_report.c)
target_include_directories(model_report PRIVATE \${PROJECT_SOURCE_DIR}/include)
target_compile_options(model_report PRIVATE${FEAT})
target_link_libraries(model_report PRIVATE solver_double PkgConfig::RAYLIB \${LIBM} \${LIBPTH} \${LIBDL} \${LIBGL})
CMAKE
        echo "features:${FEAT:- none (baseline tree)}"
    fi
fi

REV="$(git rev-parse --short "$REF")"
echo "$REV" > "$DEST/.sandbox_commit"
echo "sandbox:  $DEST"
echo "commit:   $REV"
echo "seeded:   ${#SEED[@]} data/config files + config/ablation"
echo
echo "build:    mkdir -p /share && ln -sfn /usr/share/cmake /share/cmake && \\"
echo "          cmake -S $REPO/$DEST -B $REPO/$DEST/bld -G Ninja -DCMAKE_BUILD_TYPE=Release && \\"
echo "          cmake --build $REPO/$DEST/bld --target model_report -j16"
echo "run:      cd $REPO/$DEST && ./bld/<tool>      # cwd must be the sandbox root"
echo "remove:   rm -rf $DEST"
