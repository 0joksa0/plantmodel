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

        o->result->sucrose[o->step] = tmp.carbohydrates.sucrose;
        o->result->starch[o->step] = tmp.carbohydrates.starch;
        o->result->ph[o->step] = tmp.photo.photosynthesis;
        o->result->partition[o->step] = tmp.growth.total_biomass;

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

    int is_light = (tmp.photo.light != REAL(0.0));

    if (ctx->was_light && !is_light) {
        ctx->starch_night_start = x[observer_ctx->starch_index];
    }

    ctx->was_light = is_light;
}
