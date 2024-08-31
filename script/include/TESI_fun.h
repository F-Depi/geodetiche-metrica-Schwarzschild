#ifndef TESI_FUN
#define TESI_FUN

#define R_MAX 1000

double TESI_fun_r(double r, double E, double l, int *sign, int *Nturns);


double TESI_fun_phi(double r, double l);


double TESI_fun_t(double r, double E);


int TESI_RK4(double h, double tau, double *r, double *phi, double *t,
                                    double E, double l, int *sign, int *Nturns);


double TESI_Veff(double r, double l);


void TESI_Veff_max_min(double l, double *r_max_min);


double TESI_bisezione(double a, double b, double l, double E);


void TESI_turning_points(double l, double E, double *r12);


// Used in check_parameters() for l < sqrt(3)
void TESI_m_case1(double l, double E, double *r0, double *r_lim, int *sign);


// Used in check_parameters() for sqrt(3) < l < 3
void TESI_m_case2(double l, double E, double *r0, double *r_lim, int *sign);

    
// Used in check_parameters() for l > 3
void TESI_m_case3(double l, double E, double *r0, double *r_lim, int *sign);


void check_parameters(double l, double E, double *r0, double *r_lim, int *sign);


#endif
