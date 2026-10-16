#ifndef UTIL_FPCLASS_H
#define UTIL_FPCLASS_H

/*
 * Fold-proof floating-point classification.
 *
 * The project compiles with `-ffast-math` (see add_compile_options in
 * CMakeLists.txt), which implies `-ffinite-math-only`. Under that flag the
 * compiler is entitled to assume no operand is ever NaN or Inf, and it
 * therefore constant-folds `isnan(v)` / `isinf(v)` to 0 and `isfinite(v)`
 * to 1 — the branch is removed at compile time. Every guard written with
 * those macros becomes dead code in a production build while continuing to
 * work in a plain -O2 build, so the failure is invisible in debugging.
 *
 * Observed consequence: a divergent trajectory (NaN in the nitrogen-affinity
 * state) passed the simulation NaN guard untouched and was then silently
 * recorded as zero growth by an `if (biomass > 0.0)` test, because NaN
 * compares false against everything.
 *
 * The test below reads the IEEE-754 exponent field through an integer copy.
 * An all-ones exponent means Inf or NaN. The arithmetic is integer, so no
 * fast-math assumption applies and the branch survives optimisation.
 */

#include <stdint.h>
#include <string.h>

static inline int pm_isfinite(double v)
{
    uint64_t bits;
    memcpy(&bits, &v, sizeof bits);
    return ((bits >> 52) & 0x7FFU) != 0x7FFU;
}

/* Convenience negation; reads better at call sites that test for failure. */
static inline int pm_isnonfinite(double v)
{
    return !pm_isfinite(v);
}

#endif /* UTIL_FPCLASS_H */
