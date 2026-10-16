#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
static int pm_isnonfinite(double x){uint64_t u;memcpy(&u,&x,8);return ((u>>52)&0x7FF)==0x7FF;}
double zero(void){return 0.0;}
int main(void){
  double n = zero()/zero();
  printf("isnan(n)=%d isfinite(n)=%d pm_isnonfinite(n)=%d val=%f\n", isnan(n), isfinite(n), pm_isnonfinite(n), n);
  return 0;
}
