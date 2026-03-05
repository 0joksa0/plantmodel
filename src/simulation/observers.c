#include "simulation/observers.h"
#include "model/model.h"

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

        o->result->sucrose[o->step] = tmp.core.sucrose;
        o->result->starch[o->step] = tmp.core.starch;
        o->result->ph[o->step] = tmp.core.photosynthesis;
        o->result->partition[o->step] = tmp.core.total_biomass;

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
        (double)tmp.core.light,
        (double)tmp.core.photosynthesis,
        (double)tmp.core.leaf_biomass,
        (double)tmp.core.root_biomass,
        (double)tmp.core.total_biomass,
        (double)tmp.core.sucrose,
        (double)tmp.core.starch,
        (double)tmp.core.nitrogen,
        (double)tmp.core.phosphorus,
        (double)tmp.core.nitrogen_uptake,
        (double)tmp.core.phosphorus_uptake,
        (double)tmp.core.uptake_cost,
        (double)tmp.core.transport_cost,
        (double)tmp.core.starch_degradation_rate,
        (double)tmp.core.night_efficiency_starch);
}

void nan_guard_observer(real_t t, const real_t* x, size_t n, void* ctx)
{
    (void)ctx;
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(x[i])) {
            fprintf(stderr,
                "NaN detected at t=%f in state[%zu]\n",
                (double)t, i);
            abort();
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

    int is_light = (tmp.core.light != REAL(0.0));

    if (ctx->was_light && !is_light) {
        ctx->starch_night_start = x[observer_ctx->starch_index];
    }

    ctx->was_light = is_light;
}
