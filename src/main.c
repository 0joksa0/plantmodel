#include <stdio.h>
#include "model/input.h"

void simulate_days(int days, Input* input); 

int main() {
    printf("=== Simulacija rasta biljke (24h) ===\n");
    Input input = generate_input();
    simulate_days(16, &input);
    return 0;

}

