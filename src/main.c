


#include <stdio.h>
#include "../include/model.h"

void simulate_days(int days); // iz simulation.c

int main() {
    printf("=== Simulacija rasta biljke (24h) ===\n");
    simulate_days(365);
    return 0;
}

