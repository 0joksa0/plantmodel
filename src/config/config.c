#include "config/config.h"

#include <ctype.h>
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct {
    const char* name;
    size_t offset;
} InputFieldMap;

static const InputFieldMap INPUT_FIELDS[] = {
#include "input_field_map.inc"
};

static char* trim(char* s)
{
    while (*s && isspace((unsigned char)*s))
        s++;

    if (*s == '\0')
        return s;

    char* end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) {
        *end = '\0';
        end--;
    }
    return s;
}

static int parse_bool(const char* value, int* out)
{
    if (!value || !out)
        return -1;

    if (strcasecmp(value, "1") == 0 || strcasecmp(value, "true") == 0 || strcasecmp(value, "yes") == 0 || strcasecmp(value, "on") == 0) {
        *out = 1;
        return 0;
    }
    if (strcasecmp(value, "0") == 0 || strcasecmp(value, "false") == 0 || strcasecmp(value, "no") == 0 || strcasecmp(value, "off") == 0) {
        *out = 0;
        return 0;
    }

    return -1;
}

static int parse_real(const char* value, real_t* out)
{
    if (!value || !out)
        return -1;

    errno = 0;
    char* end = NULL;
    double d = strtod(value, &end);
    if (errno != 0 || end == value || *trim(end) != '\0')
        return -1;

    *out = (real_t)d;
    return 0;
}

static int parse_int(const char* value, int* out)
{
    if (!value || !out)
        return -1;

    errno = 0;
    char* end = NULL;
    long v = strtol(value, &end, 10);
    if (errno != 0 || end == value || *trim(end) != '\0')
        return -1;

    *out = (int)v;
    return 0;
}

static int apply_simulation_key(const char* key, const char* value, SimulationConfig* sim)
{
    if (strcmp(key, "simulation.days") == 0)
        return parse_int(value, &sim->days);
    if (strcmp(key, "simulation.dt_hours") == 0)
        return parse_real(value, &sim->dt_hours);
    if (strcmp(key, "simulation.solver_tolerance") == 0)
        return parse_real(value, &sim->solver_tolerance);
    if (strcmp(key, "simulation.gui_enabled") == 0)
        return parse_bool(value, &sim->gui_enabled);
    if (strcmp(key, "simulation.write_csv") == 0)
        return parse_bool(value, &sim->write_csv);
    if (strcmp(key, "simulation.solver_type") == 0) {
        if (strcasecmp(value, "RK4") == 0) {
            sim->solver_type = SOLVER_RK4;
            return 0;
        }
        if (strcasecmp(value, "ODE45") == 0) {
            sim->solver_type = SOLVER_ODE45;
            return 0;
        }
        if (strcasecmp(value, "RKF78") == 0) {
            sim->solver_type = SOLVER_RKF78;
            return 0;
        }
        return -1;
    }

    return 1;
}

static int apply_input_key(const char* key, const char* value, Input* input)
{
    if (strncmp(key, "input.", 6) != 0)
        return 1;

    const char* field = key + 6;
    for (size_t i = 0; i < sizeof(INPUT_FIELDS) / sizeof(INPUT_FIELDS[0]); ++i) {
        if (strcmp(field, INPUT_FIELDS[i].name) == 0) {
            real_t parsed = REAL(0.0);
            if (parse_real(value, &parsed) != 0)
                return -1;
            real_t* dst = (real_t*)((char*)input + INPUT_FIELDS[i].offset);
            *dst = parsed;
            return 0;
        }
    }

    return 1;
}

int config_load_file(const char* path, SimulationConfig* sim, Input* input)
{
    if (!path || !sim || !input)
        return -1;

    FILE* f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Could not open config file: %s\n", path);
        return -1;
    }

    char line[1024];
    int line_no = 0;
    while (fgets(line, sizeof(line), f)) {
        line_no++;
        char* s = trim(line);
        if (*s == '\0' || *s == '#' || *s == ';')
            continue;

        char* eq = strchr(s, '=');
        if (!eq) {
            fprintf(stderr, "Config parse error at line %d: missing '='\n", line_no);
            fclose(f);
            return -1;
        }

        *eq = '\0';
        char* key = trim(s);
        char* value = trim(eq + 1);

        int rc = apply_simulation_key(key, value, sim);
        if (rc == 1)
            rc = apply_input_key(key, value, input);

        if (rc == 1) {
            fprintf(stderr, "Config warning line %d: unknown key '%s'\n", line_no, key);
            continue;
        }
        if (rc != 0) {
            fprintf(stderr, "Config parse error at line %d for key '%s'\n", line_no, key);
            fclose(f);
            return -1;
        }
    }

    fclose(f);
    return 0;
}

int config_get_input_field_offset(const char* field_name, size_t* offset_out)
{
    if (!field_name || !offset_out)
        return -1;

    const char* name = field_name;
    if (strncmp(name, "core.", 5) == 0) {
        name += 5;
    } else if (strncmp(name, "gas_exchange.", 13) == 0) {
        name += 13;
    } else if (strncmp(name, "input.", 6) == 0) {
        name += 6;
    }

    for (size_t i = 0; i < sizeof(INPUT_FIELDS) / sizeof(INPUT_FIELDS[0]); ++i) {
        if (strcmp(name, INPUT_FIELDS[i].name) == 0) {
            *offset_out = INPUT_FIELDS[i].offset;
            return 0;
        }
    }

    return -1;
}
