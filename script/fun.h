#ifndef FUN_H
#define FUN_H
#define P0 150.33046048             // = E_0, MeV/fm^3
#define R0 19.996542543             // km
#define M0 13.542058427             // solar masses
#define N0 0.16                     // fm^-3
#define A (13.4 * N0 / P0)          // energy density parameter
#define B (5.62 * N0 / P0)          // energy density parameter  
#define ALPHA 0.514                 // energy density exponent
#define BETA 2.436                  // energy density exponent

/*
 The system to solve is
 dm/dr = f_m = r^2 E(rho)
 dP/dr = f_P = - (P + E)(m + r^3 P)/(r^2 - 2mr)

 P(rho) = ALPHA A rho^(ALPHA+1) + BETA B rho^(BETA+1)
 rho(P) is found numerically
*/


// f_m = r^2 E
double fun_m(double r, double P, int tipo_politropica);


// f_P = - (P + E)(m + r^3 P)/(r^2 - 2mr)
double fun_P(double r, double P, double m, int tipo_politropica);


// This is used in findRho() to do Newton-Rapson
double P_of_rho(double rho);


// This is used in findRho() to do Newton-Rapson
// It's the analytic derivative of the function above
double DP_of_rho(double rho);


// Since E = E(rho), we need to invert numerically P(rho) to find rho given a P
double findRho(double P);


// Evaluates the energy for 3 different matter state functions
// tipo_politropica = 1, 2, 3
double fun_E(double P, int tipo_politropica);


// RK4 algorithm to advance 1 step of length h, in a system like
// dP/dr = fun_P
// dm/dr = fun_m 
// Updates the values of m, P given
// tipo_politropica = 1, 2, 3
void rungeKutta4(double h, double r, double *P, double *m, int tipo_politropica);

#endif
