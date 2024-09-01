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
                   double E, double l, int *sign, int *Nturns, double *r12);


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
void TESI_m_case1(double l, double E, double *r0, double *r_lim, int *sign, double *r12);


/*
 Only for sqrt(3) < l < 2, function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case2(double l, double E, double *r0, double *r_lim, int *sign, double *r12);

    
/*
 Only for l > 2, function used in check_parameters()
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void TESI_m_case3(double l, double E, double *r0, double *r_lim, int *sign, double *r12);


/*
 Predicts the orbit type from l, E, sign and r0
 Checks if r0 and sign are compatible with the l and E values otherwise changes
 r0 to the closest allowed value and sign to the correct value
 calls TESI_m_case1, TESI_m_case2 or TESI_m_case3 depending on l
*/
void check_parameters(double l, double E, double *r0, double *r_lim, int *sign,
                                                                double *r12);


#endif
