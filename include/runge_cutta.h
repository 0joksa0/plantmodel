typedef float (*ODEFunction)(float t, float y, void* params);

float runge_kutta_4(float y_n, float t_n, float h, ODEFunction f, void* params);
