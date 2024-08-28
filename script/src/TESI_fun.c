#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/TESI_fun.h"

double TESI_fun_r(double r, double E, double l, int *sign){
    double foo = 2 * (E - TESI_Veff(r, l));
    if (foo < 0){
        printf("Turning point reached\n");
        *sign *= -1;
    }
    return *sign * pow(fabs(foo), 1. / 2.);
}


double TESI_fun_phi(double r, double l){
    return l / (r * r);
}


double TESI_fun_t(double r, double E){
    return sqrt(2 * E + 1) * r / (r - 1.);
}


// RK4 algorithm to advance 1 step of length h, in a system like
// dr/d(tau) = fun_l(r, r, l)
// d(phi)/d(tau) = fun_phi(r, l)
// dt/d(tau) = fun_t(r, e)

int TESI_RK4(double h, double tau, double *r, double *phi, double *t,
                                            double E, double l, int *sign){

    double k_1 = h * TESI_fun_r(*r, E, l, sign);
    double l_1 = h * TESI_fun_phi(*r, l);
    double g_1 = h * TESI_fun_t(*r, E);

    double k_2 = h * TESI_fun_r(*r + k_1 / 2, E, l, sign);
    double l_2 = h * TESI_fun_phi(*r + k_1 / 2, l);
    double g_2 = h * TESI_fun_t(*r + k_1 / 2, E);

    double k_3 = h * TESI_fun_r(*r + k_2 / 2, E, l, sign);
    double l_3 = h * TESI_fun_phi(*r + k_2 / 2, l);
    double g_3 = h * TESI_fun_t(*r + k_2 / 2, E);

    double k_4 = h * TESI_fun_r(*r + k_3, E, l, sign);
    double l_4 = h * TESI_fun_phi(*r + k_3, l);
    double g_4 = h * TESI_fun_t(*r + k_3, E);

    *r += (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6.;
    *phi += (l_1 + 2 * l_2 + 2 * l_3 + l_4) / 6.;
    *t += (g_1 + 2 * g_2 + 2 * g_3 + g_4) / 6.;

    return 0;
}


double TESI_Veff(double r, double l){
    double foo = 1. / r;
    return (pow(l * foo, 2) * (1. - foo) - foo) / 4.;
}


void TESI_Veff_max_min(double l, double *r_12_max_min){
    
    if (l < sqrt(3)){
        r_12_max_min[0] = 0;
        r_12_max_min[1] = 0;
        r_12_max_min[2] = 0;
        r_12_max_min[3] = 0;
        printf("\nNo stationary points\n");
    }
    else if (l == sqrt(3)){
        r_12_max_min[0] = 3;
        r_12_max_min[1] = 3;
        r_12_max_min[2] = - 1. / 36.;
        r_12_max_min[3] = r_12_max_min[2];
        printf("\nr_ISCO\n");
    }
    else {
        r_12_max_min[0] = l * l * (1 - sqrt(1 - 3. / (l * l)));
        r_12_max_min[1] = l * l * (1 + sqrt(1 - 3. / (l * l)));
        r_12_max_min[2] = TESI_Veff(r_12_max_min[0], l); 
        r_12_max_min[3] = TESI_Veff(r_12_max_min[1], l); 
    }
}
