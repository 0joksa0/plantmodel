# Standalone numerical probes

Four single-file programs that reproduce defects found during the revision.
They are kept because each one isolates a failure mode that the regular
regression suite cannot see: the observable output is normal in every case, so
only a purpose-built check exposes the fault.

| File | Reproduces | Expected verdict |
|---|---|---|
| `integrator_consistency.c` | `SOLVER_RKF78` in the external solver library is inconsistent: 10.3 % relative error on `dy/dt = -y` independent of tolerance, implied decay rate 0.902380 = 1 - 82/840. | `RKF78 consistent = NO`; `ODE45` and `RK4` pass. Not a plantmodel defect; production uses `ODE45`. |
| `integrator_stall_threshold.c` | Bisects the `lambda_sdi` value at which the adaptive step collapses without the state going non-finite. At 12 h the threshold is ~0.5163; below it the run completes normally. | `pp=12 threshold in [0.5162, 0.5165]`; no stall at 4, 6, 8 h. Fixed by `SIMULATION_MAX_SOLVER_STEPS`. |
| `nan_guard_under_fastmath.c` | Under production flags (`-O3 -ffast-math`) `isnan(NaN)` returns 0 and `isfinite(NaN)` returns 1, so a NaN guard written with those predicates is dead code. | `-O2`: `isnan=1`. `-ffast-math`: `isnan=0` but `pm_isnonfinite=1`. Fixed by the bit-pattern test in `include/util/fpclass.h`. |
| `state_bound_audit.c` | Compares each state against its nominal bound over a full run. `max_starch` is recomputed algebraically as `6*(24-P)`, so starch never exceeds its effective cap; `max_sucrose` is not recomputed, and cytosolic sucrose exceeds 2.0 in up to 38.8 % of samples at 12 h. | starch `above=0` at every photoperiod; sucrose above cap 0.3 / 7.7 / 19.0 / 38.8 % at 4 / 6 / 8 / 12 h. |

## Building

Each file is standalone and links the production model sources:

```sh
gcc -O2 -DSOLVER_REAL_AS=double -Iinclude -I/usr/local/include \
    tests/numerical/<name>.c src/model/*.c src/config/*.c \
    src/simulation/*.c src/export/*.c \
    -o /tmp/<name> -lmodels_double -lsolver_double -lm
```

`nan_guard_under_fastmath.c` must be built twice, once with `-O2` and once with
`-O3 -ffast-math`, since the point of the test is the difference between them.
Run from the repository root: the runners resolve `config/` relatively.
