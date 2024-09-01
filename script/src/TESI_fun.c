#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/TESI_fun.h"


/****************************** System functions ******************************/

// dr / dtau = sqrt(2E - 2Veff)
// If the argument is negative, change the sign and increase Nturns by 1
double TESI_fun_r(double r, double E, double l, int *sign, int *Nturns){
    double foo = 2 * (E - TESI_Veff(r, l));
    if (foo < 0){
        *sign *= -1;
        *Nturns += 1;
        return 1;
    }
    return *sign * pow(foo, 1. / 2.);
}


// dphi / dtau = l / r^2
double TESI_fun_phi(double r, double l){
    return l / (r * r);
}


// dt / dtau = sqrt(2E + 1) * r / (r - 1)
// If r < 1, return 0
double TESI_fun_t(double r, double E){
    if (r < 1)
        return 0;
    return sqrt(2 * E + 1) * r / (r - 1.);
}


// Veff(r) = l^2 / (2r^2) - 1 / (2r) - l^2 / (2r^3)
double TESI_Veff(double r, double l){
    double foo = 1. / r;
    return ((l * l * foo * foo) * (1. - foo) - foo) / 2.;
    //return (l * l * (r - 1.) - r * r) / (2 * r * r * r);
}
/******************************************************************************/


/*
 h is the step size
 tau, r, phi, t are the system coordinates
 l, E are the system parameters
 sign is the sign of dr / dtau
 Nturns is the number of turns the particle has done

 RK4 algorithm to advance 1 step of length h, in a system like
 dr/d(tau) = fun_l(r, r, l, sign, Nturns)
 d(phi)/d(tau) = fun_phi(r, l)
 dt/d(tau) = fun_t(r, E)
*/
int TESI_RK4(double h, double tau, double *r, double *phi, double *t,
                                    double E, double l, int *sign, int *Nturns){

    int prev_sign = *sign;

    double k_1 = h * TESI_fun_r(*r, E, l, sign, Nturns);
    double l_1 = h * TESI_fun_phi(*r, l);
    double g_1 = h * TESI_fun_t(*r, E);

    if (*sign != prev_sign){
        printf("\nsign changed\n");
        return 1;
    }

    double k_2 = h * TESI_fun_r(*r + k_1 / 2, E, l, sign, Nturns);
    double l_2 = h * TESI_fun_phi(*r + k_1 / 2, l);
    double g_2 = h * TESI_fun_t(*r + k_1 / 2, E);

    if (*sign != prev_sign){
        printf("\nsign changed\n");
        return 1;
    }

    double k_3 = h * TESI_fun_r(*r + k_2 / 2, E, l, sign, Nturns);
    double l_3 = h * TESI_fun_phi(*r + k_2 / 2, l);
    double g_3 = h * TESI_fun_t(*r + k_2 / 2, E);

    if (*sign != prev_sign){
        printf("\nsign changed\n");
        return 1;
    }

    double k_4 = h * TESI_fun_r(*r + k_3, E, l, sign, Nturns);
    double l_4 = h * TESI_fun_phi(*r + k_3, l);
    double g_4 = h * TESI_fun_t(*r + k_3, E);

    if (*sign != prev_sign){
        printf("\nsign changed\n");
        return 1;
    }

    *r += (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6.;
    *phi += (l_1 + 2 * l_2 + 2 * l_3 + l_4) / 6.;
    *t += (g_1 + 2 * g_2 + 2 * g_3 + g_4) / 6.;

    return 0;
}


/*
 This function find a specific spet size h_tuned such that the particle gets
 at the turning point r12[0] +/- dR_MIN, then uses it to evolve all the
 variables. It should only be called when the particle is as close as possible
 to a turning point and TESI_fun_r in TESI_RK4 is changing the variable sign
 because it detects a negative argument in the square root.
 Give sign as changed by the previous call to TESI_RK4 so that this function
 knows if the particle needs to go inwards or outwards.
 h is the step size use in the RK4 algorithm
 tau, r, phi, t are the system coordinates
 l, E are the system parameters
 sign is the sign of dr / dtau
 Nturns is the number of turns the particle has done
 r12 is a vector of length 2 with theturning points r12[0] < r2[0]
 if there are no turning points this function should not be called
 if there is just one turning point r12[0] = r12[1]

 dr/d(tau) = fun_l(r, r, l, sign, Nturns)
 d(phi)/d(tau) = fun_phi(r, l)
 dt/d(tau) = fun_t(r, E)
*/
int TESI_dynamic_h(double h, double tau, double *r, double *phi, double *t,
                   double E, double l, int *sign, int *Nturns, double *r12){

    // Determine wich turning point is closer
    int i;
    if (fabs(*r - r12[0]) <= fabs(*r - r12[1]))
        i = 0;
    else
        i = 1;
    
    double h_tuned = (r12[i] - *r + (*sign) * dR_MIN)
                   / fabs(TESI_fun_r(*r, E, l, sign, Nturns));
    if (h_tuned > h) {
        printf("h_tuned > h, %.3e\n", h_tuned);
        exit(1);
    }
    if (h_tuned < 0) {
        printf("h_tuned < 0, %.3e\n", h_tuned);
        exit(1);
    }
    printf("h_tuned = %.3e\n", h_tuned);

    return 0;
}


/*
 Finds the points where Veff(r) has a maximum and a minimum r_max < r_min
 respectively.
 If there is just one stationary point, r_max = r_min.
 If there are no stationary points, r_max = r_min = 0.
 3rd and 4th elements of r_max_min are Veff(r_max) and Veff(r_min) respectively.
*/
void TESI_Veff_max_min(double l, double *r_max_min){
    
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


// Good old bisection method, already adapted to find zeros of E - Veff(r)
// in the interval [a, b]
double TESI_bisezione(double a, double b, double l, double E){
    double start = a;
    double end = b;
    double mid;

    if ((E - TESI_Veff(start, l)) * (E - TESI_Veff(end, l)) > 0){
        printf("No turning points in the interval\n");
        exit(1);
    }

    int kk = 0;
    while (fabs(start - end) > dR_MIN){

        mid = (start + end) / 2;
        double control = (E - TESI_Veff(start, l)) * (E - TESI_Veff(mid, l));

        if (control < 0) {
            end = mid;
        } else {
            start = mid;
        }
        kk++;
    }
    return mid;
}


/*
 Uses bisection method to find the turning points.
 A turning point is a point where dr/dt = 0 that implies E - Veff(r) = 0
 r12[0] = r1 is the inner turning point
 r12[1] = r2 is the outer turning point
*/
void TESI_turning_points(double l, double E, double *r12){

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
        printf("r2 too big (r2 > 1e3), choose a smaller E or l\n");
        exit(1);
    }
    r12[1] = TESI_bisezione(a, b, l, E);
}


/***************** User input handling and Orbit Recognition ******************/

/*
 Only for l < sqrt(3), function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 r12 is a vector of length 2 with the turning points r12[0] < r2[0],
    if there is 1 turning point r12[0] = r12[1],
    if there are no turning points r12[0] = r12[1] = 0.
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case1(double l, double E, double *r0, double *r_lim, int *sign,
                                                                double *r12){

    printf("Veff(r) has no local maxima or minima\n");
    // Positive E
    if (E >= 0){
        printf("E\t%.3f, ", E);

        if (*sign == 1){
            printf("positive solution chosen (sign = %d) => escape\n", *sign);
            *r_lim = 10 * (*r0);
        }
        else if (*sign == -1)
            printf("negative solution chosen (sign = %d) => infall\n", *sign);

        printf("r0\t%.3f\n", *r0);
    }

    // Negative E
    else {
        printf("E\t%.3f => infall\n", E);
        r12[0] = TESI_bisezione(1, R_MAX, l, E);
        r12[1] = r12[0]; // There is only one turning point
        if (*r0 >= *r_lim){
            *r0 = (1 - dR_MIN) * (*r_lim);
            *sign = -1;
            printf("r0\t%.3f\t(That's the maximum allowed)\n", *r0);
            printf("sign\t%d\t(r0 it's the extreme so the negative solution is forced)\n", *sign);
        }
        else{
            printf("r0\t%.3f\n", *r0);
            printf("sign = %d\n", *sign);
        }
    }
}


/*
 Only for sqrt(3) < l < 2, function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 r12 is a vector of length 2 with the turning points r12[0] < r2[0],
    if there is 1 turning point r12[0] = r12[1],
    if there are no turning points r12[0] = r12[1] = 0.
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case2(double l, double E, double *r0, double *r_lim, int *sign,
                                                                double *r12){

    double V_data[4];
    TESI_Veff_max_min(l, V_data);
    double r_max = V_data[0];
    double r_min = V_data[1];
    double Vmax = V_data[2];
    double Vmin = V_data[3];
    printf("V_max\t%.6e\t\tr_max\t%.6e\n", V_data[2], V_data[0]);
    printf("V_min\t%.6e\t\tr_min\t%.6e\n", V_data[3], V_data[1]);

    // Positive E
    if (E >= 0){
        printf("E\t%.3f, ", E);

        if (*sign == 1){
            printf("=> escape (positive solution chosen)\n");
            printf("sign\t%d\n", *sign);
            *r_lim = 10 * (*r0);
        }
        else if (*sign == -1){
            printf("=> infall (negative solution chosen)\n");
            printf("sign\t%d\n", *sign);
        }

        printf("r0\t%.3f\n", *r0);
    }

    // Vmax < E < 0
    else if (E > Vmax && E < 0){
        printf("E\t%.3f => infall\n", E);
        *r_lim = TESI_bisezione(1, R_MAX, l, E);
        r12[0] = *r_lim; // There is only one outer turning point
        r12[1] = *r_lim; // There is only one outer turning point
        if (*r0 >= *r_lim){
            *r0 = (1 - dR_MIN) * (*r_lim);
            *sign = -1;
            printf("r0\t%.3f\t(That's the maximum allowed)\n", *r0);
            printf("sign\t%d\t(r0 it's the extreme so the negative solution is forced)\n", *sign);
        }
        else{
            printf("r0\t%.3f\n", *r0);
            printf("sign = %d\n", *sign);
        }
    }
    // Check if a small r0 was selected with the intention to be beyond r_max
    // If the user wanted to be beyond r_max
    if (*r0 < r_max || E < Vmin){
        printf("E\t%.3f => infall", E);
        if (*r0 < r_max)
            printf(" (r0 < r_max was chosen)\n");
        *r_lim = TESI_bisezione(1, r_max, l, E);
        r12[0] = *r_lim; // There is only one outer turning point
        r12[1] = *r_lim; // There is only one outer turning point

        if (*r0 >= *r_lim){
            *r0 = (1 - dR_MIN) * *r_lim;
            *sign = -1;
            printf("r0\t%.3f\t(That's the maximum allowed)\n", *r0);
            printf("sign\t%d\t(r0 is the extreme so the negative solution is forced)\n", *sign);
        }
        else
            printf("r0\t%.3f\t\n", *r0);
    }

    // Vmin < E < Vmax (r1 < r0 < r2)
    else if (E < Vmax && E > Vmin){
        TESI_turning_points(l, E, r12);
        *r_lim = r12[1];
        printf("E\t%.3f => bound orbit, r1 = %.3f, r2 = %.3f\n", E, r12[0], r12[1]);
        if (*r0 < r12[0]){
            *r0 = (1 + dR_MIN) * r12[0];
            *sign = 1;
            printf("r0\t%.3f\t(That's the minimum allowed)\n", *r0);
            printf("sign\t%d\t(r0 = r1 so the positive solution is forced)\n", *sign);
        }
        else if (*r0 > r12[1]){
            *r0 = (1 - dR_MIN) * r12[1];
            *sign = -1;
            printf("r0\t%.3f\t(That's the maximum allowed)\n", *r0);
            printf("sign\t%d\t(r0 = r2 so the negative solution is forced)\n", *sign);
        }
        else {
            printf("r0\t%.3f\n", *r0);
            printf("sign\t%d\n", *sign);
        }
    }
    else {
        printf("I don't know this case from TESI_m_case2\n");
        exit(1);
    }
}


/*
 Only for l > 2, function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 r12 is a vector of length 2 with the turning points r12[0] < r2[0],
    if there is 1 turning point r12[0] = r12[1],
    if there are no turning points r12[0] = r12[1] = 0.
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case3(double l, double E, double *r0, double *r_lim, int *sign,
                                                                double *r12){

    double V_data[4];
    TESI_Veff_max_min(l, V_data);
    double r_max = V_data[0];
    double r_min = V_data[1];
    double Vmax = V_data[2];
    double Vmin = V_data[3];
    printf("V_max\t%.6e\t\tr_max\t%.6e\n", V_data[2], V_data[0]);
    printf("V_min\t%.6e\t\tr_min\t%.6e\n", V_data[3], V_data[1]);

    // E > Vmax
    if (E >= Vmax){
        printf("E\t%.3f, ", E);

        if (*sign == 1){
            printf("=> escape (positive solution chosen)\n");
            printf("sign\t%d\n", *sign);
            *r_lim = 10 * (*r0);
        }
        else if (*sign == -1){
            printf("=> infall (negative solution chosen)\n");
            printf("sign\t%d\n", *sign);
        }

        printf("r0\t%.3f\n", *r0);
    }
    
    // Check if E is so low that the particle is already beyond r_max
    // or r0 is smaller than r_max (=> the user wants to be beyond r_max)
    else if (E < Vmin || *r0 < r_max){
        printf("E\t%.3f => infall", E);
        if (*r0 < r_max)
            printf(" (r0 < r_max was chosen)\n");
        *r_lim = TESI_bisezione(1, r_max, l, E);
        r12[0] = *r_lim; // There is only one outer turning point
        r12[1] = *r_lim; // There is only one outer turning point

        if (*r0 >= *r_lim){
            *r0 = (1 - dR_MIN) * *r_lim;
            *sign = -1;
            printf("r0\t%.3f\t(That's the maximum allowed)\n", *r0);
            printf("sign\t%d\t(r0 is the extreme so the negative solution is forced)\n", *sign);
        }
        else{
            printf("r0\t%.3f\t\n", *r0);
            printf("sign\t%d\n", *sign);
        }
    }

    // 0 < E < Vmax (unbound orbit)
    else if (E < Vmax && E >= 0){
        printf("E\t%.3f => unbound orbit\n", E);
        r12[0] = TESI_bisezione(r_max, r_min, l, E);
        r12[1] = r12[0]; // There is only one inner turning point
        if (*r0 <= r12[0]){
            *r0 = (1 - dR_MIN) * r12[0];
            *sign = 1;
            printf("r0\t%.3f\t(That's the minimum allowed)\n", *r0);
            printf("sign\t%d\t(r0 it's the extreme so the positive solution is forced)\n", *sign);
        }
        else{
            printf("r0\t%.3f\n", *r0);
            printf("sign = %d\n", *sign);
        }
        if (*sign == 1)
            *r_lim = 10 * (*r0);
    }

    // Vmin < E < 0
    else if (E < 0 && E > Vmin){
        TESI_turning_points(l, E, r12);
        *r_lim = r12[1];
        printf("E\t%.3f => bound orbit, r1 = %.3f, r2 = %.3f\n", E, r12[0], r12[1]);
        if (*r0 < r12[0]){
            *r0 = (1 + dR_MIN) * r12[0];
            *sign = 1;
            printf("r0\t%.3f\t(That's the minimum allowed)\n", *r0);
            printf("sign\t%d\t(r0 = r1 so the positive solution is forced)\n", *sign);
        }
        else if (*r0 > r12[1]){
            *r0 = (1 - dR_MIN) * r12[1];
            *sign = -1;
            printf("r0\t%.3f\t(That's the maximum allowed)\n", *r0);
            printf("sign\t%d\t(r0 = r2 so the negative solution is forced)\n", *sign);
        }
        else {
            printf("r0\t%.3f\n", *r0);
            printf("sign\t%d\n", *sign);
        }
    }
    else if (E == Vmin){
        printf("E\t%.3f => circular orbit, r1 = r2 = %.3f\n", E, r_min);
        r12[0] = r_min; // Circular orbit, the 2 turning points coincide
        r12[1] = r_min; // Circular orbit, the 2 turning points coincide
        if (*r0 != r_min){
            *r0 = r_min;
            *r_lim = r_min;
            printf("r0\t%.3f\t(r0 != r1 = r2, r0 must be r1 = r2)\n", *r0);
        }
        else
            printf("r0\t%.3f\n", *r0);
    }
    else {
        printf("I don't know this case from TESI_m_case3\n");
        exit(1);
    }

}


/*
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 r12 is a vector of length 2 with the turning points r12[0] < r2[0],
    if there is 1 turning point r12[0] = r12[1],
    if there are no turning points r12[0] = r12[1] = 0.
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void check_parameters(double l, double E, double *r0, double *r_lim, int *sign,
        double *r12){

    printf("\nl\t%.3f\n", l);

    if (l >= 0      && l <= sqrt(3))
        TESI_m_case1(l, E, r0, r_lim, sign, r12);

    else if (l > sqrt(3) && l <= 2)
        TESI_m_case2(l, E, r0, r_lim, sign, r12);

    else if (l > 2)
        TESI_m_case3(l, E, r0, r_lim, sign, r12);

    else {
        printf("l must be positive\n");
        exit(1);
    }
}
