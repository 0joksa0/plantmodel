#include <solver.h>

#include "model/input.h"
#include <stdio.h>

#define PHOTO_PERIOD REAL(8.0)
void simulate_days(int days, Input* input);

int main()
{
    printf("=== Simulacija rasta biljke (24h) ===\n");
    Input input = generate_input();
    input.photoperiod = PHOTO_PERIOD;
    input.lambda_sni = REAL(0.16);
    input.lambda_sb = REAL(0.00344);        
    input.feedback_on_photosynthesis = REAL(0.82);     
    simulate_days(2, &input);
    return 0;
}
