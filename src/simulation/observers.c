#include "simulation/observers.h"
#include "model/model.h"
#include "model/state.h"   /* PlantStateIndex: X_N, X_P for the domain guard */
#include "util/fpclass.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

void observer_chain_cb(real_t t, const real_t* x, size_t n, void* vctx)
{
    ObserverChain* ch = (ObserverChain*)vctx;
    for (size_t i = 0; i < ch->count; ++i) {
        if (ch->obs[i]) {
            ch->obs[i](t, x, n, ch->ctx ? ch->ctx[i] : NULL);
        }
    }
}

void sim_sampling_observer(real_t t, const real_t* x, size_t n, void* vctx)
{
    (void)n;
    SimObsCtx* o = (SimObsCtx*)vctx;

    if (t + REAL(1e-15) < o->next_t)
        return;
    if (o->step >= o->max_steps)
        return;

    while (o->step < o->max_steps && o->next_t <= t + REAL(1e-15)) {
        Input tmp = *(o->rhs_ctx->base);
        o->unpack_state(&tmp, x);

        o->compute_algebraic(&tmp, fmod(o->next_t, REAL(24.0)), o->rhs_ctx->starch_night_start);

        o->result->sucrose[o->step]         = tmp.carbohydrates.sucrose;
        o->result->starch[o->step]          = tmp.carbohydrates.starch;
        o->result->ph[o->step]              = tmp.photo.photosynthesis;
        o->result->total_biomass[o->step]       = tmp.growth.total_biomass;
        o->result->sucrose_vacuole[o->step] = tmp.carbohydrates.sucrose_vacuole;
        o->result->sucrose_total[o->step]   = tmp.carbohydrates.sucrose
                                            + tmp.carbohydrates.sucrose_vacuole;

        o->step++;
        o->result->filled_steps = o->step;
        o->next_t += o->dt_sample;
    }
}

void csv_observer(real_t t, const real_t* x, size_t n, void* ctx)
{
    (void)n;
    CSVObserverCtx* o = (CSVObserverCtx*)ctx;
    PlantCtx* pctx = o->plant_ctx;

    Input tmp = *(pctx->base);
    o->unpack_state(&tmp, x);

    real_t hour = fmod(t, REAL(24.0));
    o->compute_algebraic(&tmp, hour, pctx->starch_night_start);

    fprintf(o->f,
        "%f,"
        "%f,"
        "%e,"
        "%e,%e,%e,"
        "%e,%e,"
        "%e,%e,"
        "%e,%e,"
        "%e,%e,"
        "%e,%e\n",
        (double)t,
        (double)tmp.photo.light,
        (double)tmp.photo.photosynthesis,
        (double)tmp.growth.leaf_biomass,
        (double)tmp.growth.root_biomass,
        (double)tmp.growth.total_biomass,
        (double)tmp.carbohydrates.sucrose,
        (double)tmp.carbohydrates.starch,
        (double)tmp.nutrients.nitrogen,
        (double)tmp.nutrients.phosphorus,
        (double)tmp.nutrients.nitrogen_uptake,
        (double)tmp.nutrients.phosphorus_uptake,
        (double)tmp.carbohydrates.uptake_cost,
        (double)tmp.carbohydrates.transport_cost,
        (double)tmp.carbohydrates.starch_degradation_rate,
        (double)tmp.carbohydrates.night_efficiency_starch);
}

/*
 * Records excursions of the nutrient pools out of the admissible quadrant.
 *
 * Purely diagnostic: it observes, never modifies the state, and never stops
 * the run.  It exists because "finite" and "meaningful" are not the same
 * condition for this model, and until now only the first was being checked.
 * See DomainGuardCtx in include/simulation/observers.h for why the
 * distinction matters to the callers that classify parameter sets.
 */
void domain_guard_observer(real_t t, const real_t* x, size_t n, void* ctx)
{
    DomainGuardCtx* g = (DomainGuardCtx*)ctx;
    if (!g)
        return;

    static const size_t guarded[] = { X_N, X_P };
    int any = 0;

    for (size_t k = 0; k < sizeof(guarded) / sizeof(guarded[0]); ++k) {
        const size_t i = guarded[k];
        if (i >= n)
            continue;
        const double v = (double)x[i];
        if (v < 0.0) {
            any = 1;
            if (!g->tripped) {
                g->tripped     = 1;
                g->t_first     = (double)t;
                g->index_first = i;
            }
            if (v < g->min_value) {
                g->min_value = v;
                g->index_min = i;
            }
        }
    }
    if (any)
        g->n_events++;
}

void nan_guard_observer(real_t t, const real_t* x, size_t n, void* ctx)
{
    NanGuardCtx* g = (NanGuardCtx*)ctx;

    for (size_t i = 0; i < n; ++i) {
        /* pm_isfinite, not isfinite: the project builds with -ffast-math,
         * under which the library macro folds to a constant 1 and this
         * entire guard is optimised away. See include/util/fpclass.h. */
        if (pm_isnonfinite((double)x[i])) {
            if (!g) {
                /* No record attached: preserve the historical fail-fast. */
                fprintf(stderr, "non-finite state at t=%f in state[%zu]\n",
                        (double)t, i);
                abort();
            }
            g->n_events++;
            if (!g->tripped) {
                g->tripped     = 1;
                g->t_first     = (double)t;
                g->index_first = i;
            }
            if (!g->warned) {
                g->warned = 1;
                fprintf(stderr,
                    "warning: trajectory diverged — non-finite state[%zu] "
                    "at t=%.3f h; integration continues, caller must classify "
                    "this run as invalid\n", i, (double)t);
            }
            return;  /* one report per sample is enough */
        }
    }
}

void night_start_observer(real_t t, const real_t* x, size_t n, void* vctx)
{
    (void)n;
    NightStartObserverCtx* observer_ctx = (NightStartObserverCtx*)vctx;
    PlantCtx* ctx = observer_ctx->plant_ctx;

    real_t hour = fmod(t, REAL(24.0));

    Input tmp = *(ctx->base);
    update_light_conditions(&tmp, hour);

    int is_light = (tmp.photo.light != REAL(0.0));

    if (ctx->was_light && !is_light) {
        ctx->starch_night_start = x[observer_ctx->starch_index];
    }

    ctx->was_light = is_light;
}
