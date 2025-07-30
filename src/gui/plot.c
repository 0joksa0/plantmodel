#include "gui/plot.h"
#include "raylib.h"
#include "raymath.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    const char* name;
    real_t* values;
    bool active;
    Color color;
    Rectangle buttonBounds;

} DataPlotEntry;

Rectangle ShrinkRect(Rectangle r, real_t amount)
{
    r.x += amount;
    r.y += amount;
    r.width -= 2 * amount;
    r.height -= 2 * amount;
    return r;
}

real_t ComputeYScale(DataPlotEntry* entries, int count, int total_steps, real_t screenHeight, real_t topPadding, real_t bottomPadding)
{
    real_t max_val = REAL(0.0001);
    for (int i = 0; i < count; i++) {
        if (!entries[i].active)
            continue;

        for (int j = 0; j < total_steps; j++) {
            if (entries[i].values[j] > max_val) {
                max_val = entries[i].values[j];
            }
        }
    }

    real_t plotHeight = screenHeight - topPadding;
    return plotHeight / max_val;
}

void DrawDataPlot(DataPlotEntry* entries, int count, int total_steps, real_t y_scale, real_t centerY, int margin)
{
    Vector2 mouse = GetMousePosition();

    for (int i = 0; i < count; i++) {
        DataPlotEntry* e = &entries[i];

        DrawRectangleRec(e->buttonBounds, e->active ? e->color : LIGHTGRAY);
        DrawText(e->name, e->buttonBounds.x + e->buttonBounds.width * REAL(1.5), e->buttonBounds.y, 20, BLACK);

        if (CheckCollisionPointRec(mouse, e->buttonBounds) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            e->active = !(e->active);
        }

        if (!e->active)
            continue;

        for (int j = 1; j < total_steps; ++j) {
            int x1 = margin + (j - 1) * (GetScreenWidth() - margin * 2) / total_steps;
            int y1 = centerY - e->values[j - 1] * y_scale;
            int x2 = margin + j * (GetScreenWidth() - margin * 2) / total_steps;
            int y2 = centerY - e->values[j] * y_scale;
            DrawLineEx((Vector2) { x1, y1 }, (Vector2) { x2, y2 }, REAL(2.0), e->color);
        }
    }
}

void printGrid(int width, int height, int margine, int steps, int total_steps, real_t y_scale)
{

    real_t tickStep = (((real_t)height / 2) - 40) / y_scale;
    tickStep /= steps;
    for (real_t i = 1; i <= steps; i++) {
        real_t val = i * tickStep;

        real_t y = ((real_t)height / 2) - (val * y_scale);
        DrawLine(margine, y, width, y, LIGHTGRAY);
        DrawText(TextFormat("%.1f", val), margine - 5 - MeasureText(TextFormat("%.1f", val), 12), y - 8, 12, DARKGRAY);
        y = ((real_t)height / 2) + (val * y_scale);

        DrawLine(margine, y, width, y, LIGHTGRAY);
        DrawText(TextFormat("-%.1f", val), margine - 5 - MeasureText(TextFormat("-%.1f", val), 12), y - 8, 12, DARKGRAY);
    }

    for (real_t i = 0; i <= total_steps; i += 24) {
        DrawLine(margine + i, 0, margine + i, height, LIGHTGRAY);
        DrawText(TextFormat("%.0f", i), margine + i - MeasureText(TextFormat("%.0f", i), 12) / 2, height / 2 + MeasureTextEx(GetFontDefault(), TextFormat("%.0f", i), 12, 1).y / 2, 12, DARKGRAY);
    }
    DrawLine(0, height / 2, width, height / 2, BLACK);
    DrawLine(margine, 0, margine, height, BLACK);
}

void main_thread(void)
{
    InitWindow(GetMonitorWidth(GetCurrentMonitor()),
        GetMonitorHeight(GetCurrentMonitor()),
        "Plant Simulation");
    SetTargetFPS(60);
    real_t screenHeight = GetMonitorHeight(GetCurrentMonitor());
    real_t screenWidht = GetMonitorWidth(GetCurrentMonitor());
    int margine = 30;
    real_t y_scale = 1;
    real_t max_val = 0;

    DataPlotEntry plots[4] = {
        { "Sucrose", sucrose, false, RED, { screenWidht - 250, 20, 20, 20 } },
        { "Starch", starch, false, BLUE, { screenWidht - 250, 50, 20, 20 } },
        { "Photos", ph, false, GREEN, { screenWidht - 250, 80, 20, 20 } },
        { "Y", partition, false, BLACK, { screenWidht - 250, 110, 20, 20 } },
    };
    while (!WindowShouldClose()) {
        /* –– PAN –– */
        SetTraceLogLevel(LOG_ERROR);

        BeginDrawing();
        ClearBackground(RAYWHITE);
        y_scale = ComputeYScale(plots, 4, total_steps, screenHeight / 2, 40, 40);
        char text[30];
        sprintf(text, "%f.0 horus", photoperiod);
        DrawText(text, screenWidht / 2, 10, 15, GREEN);
        printGrid(screenWidht, screenHeight, margine, 20, total_steps, y_scale);
        DrawDataPlot(plots, 4, total_steps, y_scale, screenHeight / 2, margine);

        EndDrawing();
    }
    CloseWindow();
}
