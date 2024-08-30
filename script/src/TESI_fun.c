#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/TESI_fun.h"

double TESI_fun_r(double r, double E, double l, int *sign, int *Nturns){
    double foo = 2 * (E - TESI_Veff(r, l));
    if (foo < 0){
        *sign *= -1;
        *Nturns += 1;
    }
    return *sign * pow(fabs(foo), 1. / 2.);
}


double TESI_fun_phi(double r, double l){
    return l / (r * r);
}


double TESI_fun_t(double r, double E){
    if (r < 1)
        return 0;
    return sqrt(2 * E + 1) * r / (r - 1.);
}


// RK4 algorithm to advance 1 step of length h, in a system like
// dr/d(tau) = fun_l(r, r, l)
// d(phi)/d(tau) = fun_phi(r, l)
// dt/d(tau) = fun_t(r, e)

int TESI_RK4(double h, double tau, double *r, double *phi, double *t,
                                    double E, double l, int *sign, int *Nturns){

    double k_1 = h * TESI_fun_r(*r, E, l, sign, Nturns);
    double l_1 = h * TESI_fun_phi(*r, l);
    double g_1 = h * TESI_fun_t(*r, E);

    double k_2 = h * TESI_fun_r(*r + k_1 / 2, E, l, sign, Nturns);
    double l_2 = h * TESI_fun_phi(*r + k_1 / 2, l);
    double g_2 = h * TESI_fun_t(*r + k_1 / 2, E);

    double k_3 = h * TESI_fun_r(*r + k_2 / 2, E, l, sign, Nturns);
    double l_3 = h * TESI_fun_phi(*r + k_2 / 2, l);
    double g_3 = h * TESI_fun_t(*r + k_2 / 2, E);

    double k_4 = h * TESI_fun_r(*r + k_3, E, l, sign, Nturns);
    double l_4 = h * TESI_fun_phi(*r + k_3, l);
    double g_4 = h * TESI_fun_t(*r + k_3, E);

    *r += (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6.;
    *phi += (l_1 + 2 * l_2 + 2 * l_3 + l_4) / 6.;
    *t += (g_1 + 2 * g_2 + 2 * g_3 + g_4) / 6.;

    return 0;
}


double TESI_Veff(double r, double l){
    double foo = 1. / r;
        return (pow(l * foo, 2) * (1. - foo) - foo) / 2.;
}


void TESI_Veff_max_min(double l, double *r_max_min){
    /*
      Finds the points where Veff(r) has a maximum and a minimum r_max < r_min
      respectively.
      If there is just one stationary point, r_max = r_min.
      If there are no stationary points, r_max = r_min = 0.
      3rd and 4th elements of r_max_min are Veff(r_max) and Veff(r_min)
    */
    
    if (l < sqrt(3)){
        r_max_min[0] = 0;
        r_max_min[1] = 0;
        r_max_min[2] = 0;
        r_max_min[3] = 0;
    }
    else if (l == sqrt(3)){
        r_max_min[0] = 3;
        r_max_min[1] = 3;
        r_max_min[2] = - 1. / 36.;
        r_max_min[3] = r_max_min[2];
    }
    else {
        r_max_min[0] = l * l * (1 - sqrt(1 - 3. / (l * l)));
        r_max_min[1] = l * l * (1 + sqrt(1 - 3. / (l * l)));
        r_max_min[2] = TESI_Veff(r_max_min[0], l); 
        r_max_min[3] = TESI_Veff(r_max_min[1], l); 
    }
}


//double TESI_outer_turning_point(double l, double E){
//    /*
//      Uses bisection method to find the turning outer turning point.
//      a turning point is a point where dr/dt = 0 that implies
//      E - Veff(r) = 0
//    */
//
//    double r_max_min[4];
//    TESI_Veff_max_min(l, r_max_min);
//
//    // We know that it must be greater than r_min
//    double a = r_max_min[1];
//    double b = 1000;
//    double c = (a + b) / 2;
//    double Vc = TESI_Veff(c, l) - E;
//    while (fabs(Vc) > 1e-8 || Vc > 0){ // E > Veff or the solver breaks
//        if (Vc < 0){
//            a = c;
//        }
//        else {
//            b = c;
//        }
//        c = (a + b) / 2;
//        Vc = TESI_Veff(c, l) - E;
//    }
//    return c;
//}
double TESI_extreme_turning_point(double l, double E){
    /*
      Uses bisection method to find the turning point for a particle already at
      r < r_max
      a turning point is a point where dr/dt = 0 that implies
      E - Veff(r) = 0
    */

    double r_max_min[4];
    TESI_Veff_max_min(l, r_max_min);

    // We know that it must be smaller than r_max
    double a = 1;
    double b = r_max_min[0];
    return TESI_bisezione(a, b, l, E);
}


void TESI_turning_points(double l, double E, double *r12){
    /*
     * Uses bisection method to find the turning points.
     * a turning point is a point where dr/dt = 0 that implies
     * E - Veff(r) = 0
     * r12[0] = r1, r12[1] = r2
     * r1 is the inner turning point
     * r2 is the outer turning point
     */

    double r_max_min[4];
    TESI_Veff_max_min(l, r_max_min);

    // We know that r_max < r1 < r_min < r2
    // Find r1
    double a = r_max_min[0];
    double b = r_max_min[1];
    r12[0] = TESI_bisezione(a, b, l, E);

    // Find r2
    a = r_max_min[1];
    b = 1000;
    if (E - TESI_Veff(b, l) > 0){
        printf("r2 too big (r2 > 1e3), choose a smoller E or l\n");
        exit(1);
    }
    r12[1] = TESI_bisezione(a, b, l, E);
}


double TESI_bisezione(double a, double b, double l, double E){
    double start = a;
    double end = b;
    double mid;

    if ((E - TESI_Veff(start, l)) * (E - TESI_Veff(end, l)) > 0){
        printf("No turning points in the interval\n");
        exit(1);
    }

    int kk = 0;
    while (fabs(start - end) > 1e-8 && kk < 100){

        mid = (start + end) / 2;
        double control = (E - TESI_Veff(start, l)) * (E - TESI_Veff(mid, l));

        if (control < 0) {
            end = mid;
        } else {
            start = mid;
        }
    }

    return mid;
}

