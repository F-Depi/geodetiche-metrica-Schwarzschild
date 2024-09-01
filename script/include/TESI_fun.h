#ifndef TESI_FUN
#define TESI_FUN

#define dR_MIN 1e-8
#define R_MAX 1000

/****************************** System functions ******************************/

// dr / dtau = sqrt(2E - 2Veff)
// If the argument is negative, change the sign and increase Nturns by 1
double TESI_fun_r(double r, double E, double l, int *sign, int *Nturns);


// dphi / dtau = l / r^2
double TESI_fun_phi(double r, double l);


// dt / dtau = sqrt(2E + 1) * r / (r - 1)
// If r < 1, return 0
double TESI_fun_t(double r, double E);


// Veff(r) = l^2 / (2r^2) - 1 / (2r) - l^2 / (2r^3)
double TESI_Veff(double r, double l);
/******************************************************************************/


// RK4 algorithm to advance 1 step of length h, in a system like
// dr/d(tau) = fun_l(r, r, l, sign, Nturns)
// d(phi)/d(tau) = fun_phi(r, l)
// dt/d(tau) = fun_t(r, E)
int TESI_RK4(double h, double tau, double *r, double *phi, double *t,
                                    double E, double l, int *sign, int *Nturns);


// Good old bisection method, already adapted to find zeros of E - Veff(r)
// in the interval [a, b]
double TESI_bisezione(double a, double b, double l, double E);


/*
 Finds the points where Veff(r) has a maximum and a minimum r_max < r_min
 respectively.
 If there is just one stationary point, r_max = r_min.
 If there are no stationary points, r_max = r_min = 0.
 3rd and 4th elements of r_max_min are Veff(r_max) and Veff(r_min) respectively.
*/
void TESI_Veff_max_min(double l, double *r_max_min);


/*
 Uses bisection method to find the turning points.
 A turning point is a point where dr/dt = 0 that implies E - Veff(r) = 0
 r12[0] = r1 is the inner turning point
 r12[1] = r2 is the outer turning point
*/
void TESI_turning_points(double l, double E, double *r12);


/***************** User input handling and Orbit Recognition ******************/

/*
 Only for l < sqrt(3), function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case1(double l, double E, double *r0, double *r_lim, int *sign);


/*
 Only for sqrt(3) < l < 2, function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case2(double l, double E, double *r0, double *r_lim, int *sign);

    
/*
 Only for l > 2, function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case3(double l, double E, double *r0, double *r_lim, int *sign);


/*
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void check_parameters(double l, double E, double *r0, double *r_lim, int *sign);


#endif
