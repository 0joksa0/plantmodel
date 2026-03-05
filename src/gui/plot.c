#include "gui/plot.h"

#include "raylib.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SERIES_COUNT 4

typedef enum {
    VIEW_FULL = 0,
    VIEW_LAST_24H,
    VIEW_LAST_48H,
    VIEW_CUSTOM
} ViewMode;

typedef struct {
    const char* name;
    const real_t* values;
    Color color;
    bool visible;
    Rectangle toggle;
} Series;

typedef struct {
    const char* label;
    ViewMode mode;
    Rectangle rect;
} ViewModeButton;

typedef struct {
    Color bg;
    Color panel;
    Color border;
    Color text;
    Color text_dim;
    Color grid_major;
    Color grid_minor;
    Color accent;
    Color day_band;
    Color night_band;
} UiTheme;

typedef struct {
    int available;
    long steps;
    long accepted;
    long rejected;
    double accept_ratio;
    double avg_dt;
    double max_err;
    double wall_time;
} SummaryStats;

static int clamp_i(int v, int lo, int hi)
{
    if (v < lo)
        return lo;
    if (v > hi)
        return hi;
    return v;
}

static real_t clamp_r(real_t v, real_t lo, real_t hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

static int valid_steps(const GuiOutput* out)
{
    if (!out || !out->filled_steps)
        return 0;
    return clamp_i(*out->filled_steps, 0, out->total_steps);
}

static int to_index(float x, Rectangle area, int from_idx, int to_idx)
{
    int span = to_idx - from_idx;
    if (span <= 0)
        return from_idx;

    float t = (x - area.x) / area.width;
    if (t < 0.0f)
        t = 0.0f;
    if (t > 1.0f)
        t = 1.0f;

    return from_idx + (int)roundf(t * (float)span);
}

static float to_x(int idx, Rectangle area, int from_idx, int to_idx)
{
    int span = to_idx - from_idx;
    if (span <= 0)
        return area.x;
    float t = (float)(idx - from_idx) / (float)span;
    return area.x + t * area.width;
}

static real_t visible_max(const Series* series, int from_idx, int to_idx)
{
    real_t maxv = REAL(0.001);
    for (int i = 0; i < SERIES_COUNT; ++i) {
        if (!series[i].visible || !series[i].values)
            continue;

        for (int j = from_idx; j <= to_idx; ++j) {
            real_t v = series[i].values[j];
            if (v > maxv)
                maxv = v;
        }
    }
    return maxv;
}

static bool any_visible(const Series* series)
{
    for (int i = 0; i < SERIES_COUNT; ++i) {
        if (series[i].visible)
            return true;
    }
    return false;
}

static Font load_ui_font(bool* custom_loaded)
{
    char path_buf[512];
    const char* candidates[] = {
        "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
        "/usr/share/fonts/truetype/noto/NotoSansMono-Regular.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
        "/usr/share/fonts/TTF/DejaVuSansMono.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf"
    };

    *custom_loaded = false;

    FILE* pipe = popen("fc-match -f '%{file}' sans 2>/dev/null", "r");
    if (pipe) {
        if (fgets(path_buf, sizeof(path_buf), pipe)) {
            path_buf[strcspn(path_buf, "\r\n")] = '\0';
            if (path_buf[0] != '\0' && FileExists(path_buf)) {
                Font f = LoadFontEx(path_buf, 64, NULL, 0);
                if (f.texture.id > 0) {
                    SetTextureFilter(f.texture, TEXTURE_FILTER_BILINEAR);
                    *custom_loaded = true;
                    pclose(pipe);
                    return f;
                }
            }
        }
        pclose(pipe);
    }

    pipe = popen("fc-match -f '%{file}' monospace 2>/dev/null", "r");
    if (pipe) {
        if (fgets(path_buf, sizeof(path_buf), pipe)) {
            path_buf[strcspn(path_buf, "\r\n")] = '\0';
            if (path_buf[0] != '\0' && FileExists(path_buf)) {
                Font f = LoadFontEx(path_buf, 64, NULL, 0);
                if (f.texture.id > 0) {
                    SetTextureFilter(f.texture, TEXTURE_FILTER_BILINEAR);
                    *custom_loaded = true;
                    pclose(pipe);
                    return f;
                }
            }
        }
        pclose(pipe);
    }

    *custom_loaded = false;
    for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
        if (!FileExists(candidates[i]))
            continue;

        Font f = LoadFontEx(candidates[i], 64, NULL, 0);
        if (f.texture.id > 0) {
            SetTextureFilter(f.texture, TEXTURE_FILTER_BILINEAR);
            *custom_loaded = true;
            return f;
        }
    }

    Font fallback = GetFontDefault();
    SetTextureFilter(fallback.texture, TEXTURE_FILTER_BILINEAR);
    return fallback;
}

static void draw_label(Font font, const char* text, Vector2 pos, float size, Color color)
{
    DrawTextEx(font, text, pos, size, 0.0f, color);
}

static int read_summary_csv(const char* path, SummaryStats* out)
{
    if (!path || !out)
        return 0;

    FILE* f = fopen(path, "r");
    if (!f)
        return 0;

    char header[512];
    char line[512];
    if (!fgets(header, sizeof(header), f) || !fgets(line, sizeof(line), f)) {
        fclose(f);
        return 0;
    }
    fclose(f);

    char* nl = strchr(line, '\n');
    if (nl)
        *nl = '\0';

    char* save = NULL;
    char* tok = strtok_r(line, ",", &save);
    int idx = 0;
    double vals[7] = { 0 };
    while (tok && idx < 7) {
        vals[idx++] = strtod(tok, NULL);
        tok = strtok_r(NULL, ",", &save);
    }
    if (idx < 7)
        return 0;

    out->steps = (long)vals[0];
    out->accepted = (long)vals[1];
    out->rejected = (long)vals[2];
    out->accept_ratio = vals[3];
    out->avg_dt = vals[4];
    out->max_err = vals[5];
    out->wall_time = vals[6];
    out->available = 1;
    return 1;
}

static void draw_summary_panel(Font font, const UiTheme* theme, Rectangle panel, const SummaryStats* s)
{
    DrawRectangleRec(panel, theme->panel);
    DrawRectangleLinesEx(panel, 1.0f, theme->border);

    draw_label(font, "solver summary", (Vector2){ panel.x + 10, panel.y + 10 }, 20.0f, theme->text);
    if (!s || !s->available) {
        draw_label(font, "summary.csv not ready", (Vector2){ panel.x + 10, panel.y + 40 }, 18.0f, theme->text_dim);
        return;
    }

    draw_label(font, TextFormat("steps: %ld", s->steps), (Vector2){ panel.x + 10, panel.y + 40 }, 18.0f, theme->text);
    draw_label(font, TextFormat("accepted: %ld", s->accepted), (Vector2){ panel.x + 10, panel.y + 64 }, 18.0f, theme->text);
    draw_label(font, TextFormat("rejected: %ld", s->rejected), (Vector2){ panel.x + 10, panel.y + 88 }, 18.0f, theme->text);
    draw_label(font, TextFormat("accept ratio: %.4f", s->accept_ratio), (Vector2){ panel.x + 10, panel.y + 112 }, 18.0f, theme->text);
    draw_label(font, TextFormat("avg dt: %.6f", s->avg_dt), (Vector2){ panel.x + 10, panel.y + 136 }, 18.0f, theme->text);
    draw_label(font, TextFormat("max err: %.3e", s->max_err), (Vector2){ panel.x + 10, panel.y + 160 }, 18.0f, theme->text);
    draw_label(font, TextFormat("wall time [s]: %.4f", s->wall_time), (Vector2){ panel.x + 10, panel.y + 184 }, 18.0f, theme->text);
}

static void draw_legend(Font font, const UiTheme* theme, Series* series)
{
    Vector2 m = GetMousePosition();

    for (int i = 0; i < SERIES_COUNT; ++i) {
        Color bg = series[i].visible ? Fade(series[i].color, 0.18f) : theme->panel;
        DrawRectangleRec(series[i].toggle, bg);
        DrawRectangleLinesEx(series[i].toggle, 1.0f, theme->border);

        DrawRectangle((int)series[i].toggle.x + 8, (int)series[i].toggle.y + 7, 12, 12, series[i].color);
        draw_label(font, series[i].name, (Vector2){ series[i].toggle.x + 28, series[i].toggle.y + 5 }, 20.0f, theme->text);

        if (CheckCollisionPointRec(m, series[i].toggle) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            series[i].visible = !series[i].visible;
        }
    }
}

static void draw_view_mode_buttons(Font font, const UiTheme* theme, ViewModeButton* buttons, ViewMode* mode)
{
    Vector2 m = GetMousePosition();

    for (int i = 0; i < 3; ++i) {
        bool active = (*mode == buttons[i].mode);
        Color bg = active ? theme->accent : theme->panel;
        Color fg = active ? RAYWHITE : theme->text;

        DrawRectangleRec(buttons[i].rect, bg);
        DrawRectangleLinesEx(buttons[i].rect, 1.0f, theme->border);
        draw_label(font, buttons[i].label, (Vector2){ buttons[i].rect.x + 14, buttons[i].rect.y + 5 }, 20.0f, fg);

        if (CheckCollisionPointRec(m, buttons[i].rect) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            *mode = buttons[i].mode;
        }
    }
}

static void draw_grid(const UiTheme* theme, Rectangle area, real_t y_min, real_t y_max, real_t t0, real_t t1)
{
    const int y_ticks = 8;
    const int x_ticks = 10;

    for (int i = 0; i <= y_ticks; ++i) {
        float k = (float)i / (float)y_ticks;
        float y = area.y + area.height - area.height * k;
        DrawLine((int)area.x, (int)y, (int)(area.x + area.width), (int)y, i == 0 ? theme->grid_major : theme->grid_minor);
        real_t y_value = y_min + ((y_max - y_min) * (real_t)k);
        DrawText(TextFormat("%.3f", (double)y_value), (int)area.x - 74, (int)y - 8, 14, theme->text_dim);
    }

    for (int i = 0; i <= x_ticks; ++i) {
        float k = (float)i / (float)x_ticks;
        float x = area.x + area.width * k;
        DrawLine((int)x, (int)area.y, (int)x, (int)(area.y + area.height), i == 0 ? theme->grid_major : theme->grid_minor);
        real_t t = t0 + (t1 - t0) * (real_t)k;
        DrawText(TextFormat("%.1f", (double)t), (int)x - 12, (int)(area.y + area.height + 6), 12, theme->text_dim);
    }

    DrawText("time [h]", (int)(area.x + area.width / 2) - 28, (int)(area.y + area.height + 24), 14, theme->text_dim);
    DrawRectangleLinesEx(area, 1.2f, theme->border);
}

static void draw_day_night(const UiTheme* theme, const GuiOutput* out, Rectangle area, int from_idx, int to_idx)
{
    if (!out || out->dt_hours <= REAL(0.0))
        return;

    int steps_per_day = (int)round((double)(REAL(24.0) / out->dt_hours));
    int light_steps = (int)round((double)(out->photoperiod / out->dt_hours));
    if (steps_per_day <= 0)
        return;

    int start_day = from_idx / steps_per_day;
    int end_day = to_idx / steps_per_day;

    for (int d = start_day; d <= end_day; ++d) {
        int day_begin = d * steps_per_day;
        int day_end = day_begin + steps_per_day;
        int light_end = day_begin + light_steps;

        int a0 = clamp_i(day_begin, from_idx, to_idx);
        int a1 = clamp_i(light_end, from_idx, to_idx);
        int a2 = clamp_i(day_end, from_idx, to_idx);

        float x0 = to_x(a0, area, from_idx, to_idx);
        float x1 = to_x(a1, area, from_idx, to_idx);
        float x2 = to_x(a2, area, from_idx, to_idx);

        if (x1 > x0)
            DrawRectangleRec((Rectangle){ x0, area.y, x1 - x0, area.height }, theme->day_band);
        if (x2 > x1)
            DrawRectangleRec((Rectangle){ x1, area.y, x2 - x1, area.height }, theme->night_band);
    }
}

static void draw_series(const Series* series, Rectangle area, int from_idx, int to_idx, real_t yscale, real_t y_min)
{
    int span = to_idx - from_idx;
    if (span <= 0)
        return;

    int draw_step = 1;
    int pixel_w = (int)area.width;
    if (pixel_w > 0 && span > pixel_w * 3) {
        draw_step = span / pixel_w;
        if (draw_step < 1)
            draw_step = 1;
    }

    for (int s = 0; s < SERIES_COUNT; ++s) {
        if (!series[s].visible || !series[s].values)
            continue;

        int i = from_idx + draw_step;
        for (; i <= to_idx; i += draw_step) {
            int prev = i - draw_step;
            float x1 = to_x(prev, area, from_idx, to_idx);
            float x2 = to_x(i, area, from_idx, to_idx);
            float y1 = area.y + area.height - (float)((series[s].values[prev] - y_min) * yscale);
            float y2 = area.y + area.height - (float)((series[s].values[i] - y_min) * yscale);
            DrawLineEx((Vector2){ x1, y1 }, (Vector2){ x2, y2 }, 2.0f, series[s].color);
        }

        if ((i - draw_step) < to_idx) {
            int prev = i - draw_step;
            float x1 = to_x(prev, area, from_idx, to_idx);
            float x2 = to_x(to_idx, area, from_idx, to_idx);
            float y1 = area.y + area.height - (float)((series[s].values[prev] - y_min) * yscale);
            float y2 = area.y + area.height - (float)((series[s].values[to_idx] - y_min) * yscale);
            DrawLineEx((Vector2){ x1, y1 }, (Vector2){ x2, y2 }, 2.0f, series[s].color);
        }
    }
}

static void draw_probe(Font font, const UiTheme* theme, const GuiOutput* out, const Series* series, Rectangle area, int from_idx, int to_idx)
{
    Vector2 m = GetMousePosition();
    if (!CheckCollisionPointRec(m, area))
        return;

    int idx = to_index(m.x, area, from_idx, to_idx);
    float x = to_x(idx, area, from_idx, to_idx);
    DrawLine((int)x, (int)area.y, (int)x, (int)(area.y + area.height), Fade(theme->text, 0.35f));

    Rectangle panel = { area.x + area.width + 16, area.y, 305, 216 };
    DrawRectangleRec(panel, theme->panel);
    DrawRectangleLinesEx(panel, 1.0f, theme->border);

    real_t hour = (real_t)idx * out->dt_hours;
    draw_label(font, TextFormat("sample %d", idx), (Vector2){ panel.x + 10, panel.y + 10 }, 20.0f, theme->text);
    draw_label(font, TextFormat("time: %.3f h", (double)hour), (Vector2){ panel.x + 10, panel.y + 38 }, 18.0f, theme->text_dim);
    draw_label(font, TextFormat("day: %d  hod: %.2f", (int)(hour / REAL(24.0)), (double)fmod(hour, REAL(24.0))), (Vector2){ panel.x + 10, panel.y + 60 }, 18.0f, theme->text_dim);

    int row = 0;
    for (int i = 0; i < SERIES_COUNT; ++i) {
        if (!series[i].visible || !series[i].values)
            continue;
        draw_label(font, TextFormat("%s: %.6f", series[i].name, (double)series[i].values[idx]), (Vector2){ panel.x + 10, panel.y + 92 + row * 24 }, 18.0f, series[i].color);
        row++;
    }
}

void* gui_main_thread(void* arg)
{
    GuiThreadArgs* args = (GuiThreadArgs*)arg;
    if (!args || !args->output)
        return NULL;

    const GuiOutput* out = args->output;
    if (!out->sucrose || !out->starch || !out->photosynthesis || !out->total_biomass || out->total_steps <= 1)
        return NULL;

#if defined(__linux__)
    if (!getenv("DISPLAY") && !getenv("WAYLAND_DISPLAY"))
        return NULL;
#endif

    SetConfigFlags(FLAG_WINDOW_HIGHDPI | FLAG_MSAA_4X_HINT);
    int width = 1520;
    int height = 900;
    InitWindow(width, height, args->title ? args->title : "PlantModel GUI");
    SetTargetFPS(args->target_fps > 0 ? args->target_fps : 60);

    bool custom_font_loaded = false;
    Font ui_font = load_ui_font(&custom_font_loaded);

    UiTheme theme = {
        .bg = { 248, 249, 251, 255 },
        .panel = { 255, 255, 255, 255 },
        .border = { 188, 194, 202, 255 },
        .text = { 24, 28, 34, 255 },
        .text_dim = { 86, 94, 105, 255 },
        .grid_major = { 176, 182, 190, 255 },
        .grid_minor = { 223, 227, 232, 255 },
        .accent = { 45, 74, 122, 255 },
        .day_band = { 242, 242, 242, 150 },
        .night_band = { 231, 235, 240, 150 }
    };

    Rectangle header = { 24, 16, (float)width - 48, 56 };
    Rectangle card = { 24, 82, 1168, 742 };
    Rectangle plot_area = { 92, 146, 1036, 602 };
    Rectangle summary_panel = { 1208, 160, 288, 232 };

    real_t zoom_y = REAL(1.0);
    real_t y_offset = REAL(0.0);
    real_t last_y_span = REAL(1.0);
    int window_hours = 48;
    ViewMode view_mode = VIEW_LAST_48H;
    int custom_center_idx = 0;
    int custom_window_steps = 0;
    int pan_active = 0;
    Vector2 pan_last_mouse = { 0 };

    Series series[SERIES_COUNT] = {
        { "Sucrose", out->sucrose, { 178, 34, 34, 255 }, true, { 92, 102, 146, 30 } },
        { "Starch", out->starch, { 49, 79, 133, 255 }, true, { 246, 102, 126, 30 } },
        { "Photosynthesis", out->photosynthesis, { 46, 115, 68, 255 }, true, { 378, 102, 196, 30 } },
        { "Total Biomass", out->total_biomass, { 30, 30, 30, 255 }, true, { 580, 102, 178, 30 } },
    };

    ViewModeButton mode_buttons[3] = {
        { "Full", VIEW_FULL, { 908, 102, 78, 30 } },
        { "24h", VIEW_LAST_24H, { 992, 102, 78, 30 } },
        { "48h", VIEW_LAST_48H, { 1076, 102, 78, 30 } },
    };
    SummaryStats summary = { 0 };
    int summary_refresh_tick = 0;

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_ONE))
            series[0].visible = !series[0].visible;
        if (IsKeyPressed(KEY_TWO))
            series[1].visible = !series[1].visible;
        if (IsKeyPressed(KEY_THREE))
            series[2].visible = !series[2].visible;
        if (IsKeyPressed(KEY_FOUR))
            series[3].visible = !series[3].visible;

        if (IsKeyPressed(KEY_A)) {
            for (int i = 0; i < SERIES_COUNT; ++i)
                series[i].visible = true;
        }
        if (IsKeyPressed(KEY_N)) {
            for (int i = 0; i < SERIES_COUNT; ++i)
                series[i].visible = false;
        }

        if (IsKeyPressed(KEY_F))
            view_mode = VIEW_FULL;
        if (IsKeyPressed(KEY_T))
            view_mode = VIEW_LAST_24H;
        if (IsKeyPressed(KEY_Y))
            view_mode = VIEW_LAST_48H;

        if (IsKeyPressed(KEY_EQUAL) || IsKeyPressed(KEY_KP_ADD)) {
            window_hours = clamp_i(window_hours + 12, 12, out->days * 24);
            view_mode = VIEW_CUSTOM;
        }
        if (IsKeyPressed(KEY_MINUS) || IsKeyPressed(KEY_KP_SUBTRACT)) {
            window_hours = clamp_i(window_hours - 12, 12, out->days * 24);
            view_mode = VIEW_CUSTOM;
        }

        if (IsKeyPressed(KEY_R))
            zoom_y = REAL(1.0);
        if (IsKeyPressed(KEY_R))
            y_offset = REAL(0.0);

        int filled = valid_steps(out);
        int full = out->total_steps;
        if (args->summary_csv_path && ((summary_refresh_tick++ % 30) == 0)) {
            read_summary_csv(args->summary_csv_path, &summary);
        }

        int from_idx = 0;
        int to_idx = filled > 1 ? filled - 1 : 0;
        int span = to_idx - from_idx;

        if (view_mode == VIEW_LAST_24H || view_mode == VIEW_LAST_48H) {
            int hours = window_hours;
            if (view_mode == VIEW_LAST_24H)
                hours = 24;
            else if (view_mode == VIEW_LAST_48H)
                hours = 48;

            int window_steps = (int)round((double)((real_t)hours / out->dt_hours));
            if (window_steps > 0 && (to_idx - from_idx) > window_steps)
                from_idx = to_idx - window_steps;
            span = to_idx - from_idx;
        } else if (view_mode == VIEW_CUSTOM && filled > 1) {
            if (custom_window_steps <= 0)
                custom_window_steps = (int)round((double)((real_t)window_hours / out->dt_hours));
            custom_window_steps = clamp_i(custom_window_steps, 8, filled - 1);
            custom_center_idx = clamp_i(custom_center_idx == 0 ? (filled - 1) : custom_center_idx, 0, filled - 1);
            from_idx = custom_center_idx - custom_window_steps / 2;
            to_idx = from_idx + custom_window_steps;
            if (from_idx < 0) {
                from_idx = 0;
                to_idx = custom_window_steps;
            }
            if (to_idx >= filled) {
                to_idx = filled - 1;
                from_idx = to_idx - custom_window_steps;
                if (from_idx < 0)
                    from_idx = 0;
            }
            span = to_idx - from_idx;
        }

        Vector2 mouse = GetMousePosition();
        int mouse_in_plot = CheckCollisionPointRec(mouse, plot_area);

        if (mouse_in_plot && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            pan_active = 1;
            pan_last_mouse = mouse;
            if (view_mode != VIEW_CUSTOM) {
                view_mode = VIEW_CUSTOM;
                custom_center_idx = (from_idx + to_idx) / 2;
                custom_window_steps = clamp_i(span, 8, filled > 1 ? filled - 1 : 8);
            }
        }
        if (pan_active && IsMouseButtonDown(MOUSE_BUTTON_LEFT) && filled > 1) {
            float dx = mouse.x - pan_last_mouse.x;
            float dy = mouse.y - pan_last_mouse.y;
            if (fabsf(dx) > 0.01f && span > 0) {
                int delta = (int)roundf((dx / plot_area.width) * (float)span);
                custom_center_idx -= delta;
                custom_center_idx = clamp_i(custom_center_idx, 0, filled - 1);
            }
            if (fabsf(dy) > 0.01f) {
                y_offset += ((real_t)dy / (real_t)plot_area.height) * last_y_span;
            }
            pan_last_mouse = mouse;
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
            pan_active = 0;

        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f && mouse_in_plot && filled > 1) {
            real_t factor = (real_t)(1.0 + wheel * 0.1);
            factor = clamp_r(factor, REAL(0.5), REAL(1.5));
            zoom_y = clamp_r(zoom_y * factor, REAL(0.2), REAL(20.0));

            if (view_mode != VIEW_CUSTOM) {
                view_mode = VIEW_CUSTOM;
                custom_center_idx = (from_idx + to_idx) / 2;
                custom_window_steps = clamp_i(span, 8, filled - 1);
            }

            int old_span = custom_window_steps > 0 ? custom_window_steps : clamp_i(span, 8, filled - 1);
            int new_span = (int)round((double)((real_t)old_span / factor));
            new_span = clamp_i(new_span, 8, filled - 1);

            int mouse_idx = to_index(mouse.x, plot_area, from_idx, to_idx);
            float rel = old_span > 0 ? ((float)(mouse_idx - from_idx) / (float)old_span) : 0.5f;
            int new_from = mouse_idx - (int)roundf(rel * (float)new_span);
            int new_to = new_from + new_span;
            if (new_from < 0) {
                new_from = 0;
                new_to = new_span;
            }
            if (new_to >= filled) {
                new_to = filled - 1;
                new_from = new_to - new_span;
                if (new_from < 0)
                    new_from = 0;
            }
            custom_window_steps = new_to - new_from;
            custom_center_idx = (new_from + new_to) / 2;
        }

        BeginDrawing();
        ClearBackground(theme.bg);

        DrawRectangleRec(header, theme.panel);
        DrawRectangleLinesEx(header, 1.0f, theme.border);
        draw_label(ui_font, "PlantModel Scientific Plot", (Vector2){ header.x + 12, header.y + 14 }, 24.0f, theme.text);
        draw_label(ui_font, "live simulation output", (Vector2){ header.x + 430, header.y + 17 }, 18.0f, theme.text_dim);

        DrawRectangleRec(card, theme.panel);
        DrawRectangleLinesEx(card, 1.0f, theme.border);

        draw_legend(ui_font, &theme, series);
        draw_view_mode_buttons(ui_font, &theme, mode_buttons, &view_mode);

        draw_label(ui_font,
            TextFormat("progress %d / %d (%.1f%%)", filled, full, full > 0 ? 100.0 * (double)filled / (double)full : 0.0),
            (Vector2){ 1210, 104 },
            18.0f,
            theme.text);
        DrawRectangle(1210, 128, 280, 10, theme.grid_minor);
        DrawRectangle(1210, 128, (int)(280.0f * (full > 0 ? (float)filled / (float)full : 0.0f)), 10, theme.accent);
        draw_summary_panel(ui_font, &theme, summary_panel, &summary);

        if (filled < 2) {
            DrawRectangleRec(plot_area, theme.panel);
            DrawRectangleLinesEx(plot_area, 1.0f, theme.border);
            draw_label(ui_font, "waiting for samples...", (Vector2){ plot_area.x + 20, plot_area.y + 24 }, 22.0f, theme.text_dim);
        } else if (!any_visible(series)) {
            DrawRectangleRec(plot_area, theme.panel);
            DrawRectangleLinesEx(plot_area, 1.0f, theme.border);
            draw_label(ui_font, "no visible series (toggle legend / keys 1-4)", (Vector2){ plot_area.x + 20, plot_area.y + 24 }, 20.0f, theme.text_dim);
        } else {
            draw_day_night(&theme, out, plot_area, from_idx, to_idx);

            real_t ymax = visible_max(series, from_idx, to_idx);
            real_t y_span = ymax / zoom_y;
            if (y_span < REAL(0.001))
                y_span = REAL(0.001);
            last_y_span = y_span;
            real_t y_min = y_offset;
            real_t y_max = y_offset + y_span;
            real_t yscale = (real_t)plot_area.height / y_span;

            real_t t0 = (real_t)from_idx * out->dt_hours;
            real_t t1 = (real_t)to_idx * out->dt_hours;
            draw_grid(&theme, plot_area, y_min, y_max, t0, t1);
            draw_series(series, plot_area, from_idx, to_idx, yscale, y_min);
            draw_probe(ui_font, &theme, out, series, plot_area, from_idx, to_idx);
        }

        const char* mode_label = (view_mode == VIEW_FULL) ? "Full" : ((view_mode == VIEW_LAST_24H) ? "24h" : ((view_mode == VIEW_LAST_48H) ? "48h" : "Custom"));
        draw_label(ui_font,
            TextFormat("photoperiod %.1fh | days %d | dt %.3fh | view %s | custom window %dh", (double)out->photoperiod, out->days, (double)out->dt_hours, mode_label, window_hours),
            (Vector2){ 32, 836 },
            18.0f,
            theme.text_dim);
        draw_label(ui_font,
            "controls: 1-4 toggle, A all, N none, F/T/Y view modes, drag to pan, wheel zoom X+Y, +/- custom window, R reset",
            (Vector2){ 32, 860 },
            16.0f,
            theme.text_dim);

        EndDrawing();
    }

    if (custom_font_loaded)
        UnloadFont(ui_font);

    CloseWindow();
    return NULL;
}
