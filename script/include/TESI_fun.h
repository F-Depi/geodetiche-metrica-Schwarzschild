#ifndef TESI_FUN
#define TESI_FUN


double TESI_fun_r(double r, double E, double l, int *sign, int *Nturns);


double TESI_fun_phi(double r, double l);


double TESI_fun_t(double r, double E);


int TESI_RK4(double h, double tau, double *r, double *phi, double *t,
                                    double E, double l, int *sign, int *Nturns);


double TESI_Veff(double r, double l);


void TESI_Veff_max_min(double l, double *r_max_min);


double TESI_bisezione(double a, double b, double l, double E);


void TESI_turning_points(double l, double E, double *r12);


double TESI_extreme_turning_point(double l, double E);


#endif
