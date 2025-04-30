#include <raylib.h>

void DrawAxisGrid(int scrollX, int scrollY, float zoom)
{
    int screenWidth = GetScreenWidth();
    int screenHeight = GetScreenHeight();

    // X-axis
    DrawLine(0, 500 - scrollY, screenWidth, 500 - scrollY, LIGHTGRAY);
    for (int i = -scrollX; i < screenWidth / zoom + scrollX; i += 24) {
        int worldX = i;
        int screenX = 20 + worldX * zoom - scrollX;
        DrawLine(screenX, 495 - scrollY, screenX, 505 - scrollY, GRAY);
        DrawText(TextFormat("%d", worldX), screenX - 10, 510 - scrollY, 10, DARKGRAY);
    }

    // Y-axis
    DrawLine(20 - scrollX, 0, 20 - scrollX, screenHeight, LIGHTGRAY);
    for (int j = -scrollY; j < screenHeight / (20 * zoom) + scrollY; j += 5) {
        int worldY = j;
        int screenY = 500 - worldY * 20 * zoom - scrollY;
        DrawLine(15 - scrollX, screenY, 25 - scrollX, screenY, GRAY);
        DrawText(TextFormat("%d", worldY), 0, screenY - 5, 10, DARKGRAY);
    }
}



void main_thread(){

      InitWindow(GetMonitorWidth(GetCurrentMonitor()), GetMonitorHeight(GetCurrentMonitor()), "Plant Simulation");
    SetTargetFPS(60);
    int scrollX = 0;
    int scrollY = 0;
    int zoom = 5.0f;
    // while (!WindowShouldClose()) {
    //     if (IsKeyDown(KEY_RIGHT))
    //         scrollX += 5;
    //     if (IsKeyDown(KEY_LEFT))
    //         scrollX -= 5;
    //     if (IsKeyDown(KEY_DOWN))
    //         scrollY += 5;
    //     if (IsKeyDown(KEY_UP))
    //         scrollY -= 5;
    //     if (IsKeyPressed(KEY_A))
    //         zoom *= 2.0f;
    //     if (IsKeyPressed(KEY_B))
    //         zoom /= 1.1f;
    //
    //     BeginDrawing();
    //     ClearBackground(RAYWHITE);
    //     DrawAxisGrid(scrollX, scrollY, zoom);
    //     for (int i = 1; i < total_steps; ++i) {
    //         DrawLine(2 + (i - 1 - scrollX) * zoom, 500 - sucroseValues[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 500 - sucroseValues[i] * 20 * zoom - scrollY, RED);
    //         DrawLine(2 + (i - 1 - scrollX) * zoom, 500 - starchValues[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 500 - starchValues[i] * 20 * zoom - scrollY, BLUE);
    //         DrawLine(2 + (i - 1 - scrollX) * zoom, 500 - phValues[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 500 - phValues[i] * 20 * zoom - scrollY, GREEN);
    //         DrawLine(2 + (i - 1 - scrollX) * zoom, 50 - partition[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 50 - partition[i] * 20 * zoom - scrollY, BLACK);
    //     }
    //
    //     DrawText("Sucrose", 20, 20, 20, RED);
    //     DrawText("Starch", 20, 40, 20, BLUE);
    //     DrawText("Photosynthesis", 20, 60, 20, GREEN);
    //     DrawText("Biomass", 20, 80, 20, BLACK);
    //     DrawText(TextFormat("Zoom: %.2fx", zoom), 1000, 20, 20, DARKGRAY);
    //
    //     EndDrawing();
    // }
    // CloseWindow();



}
