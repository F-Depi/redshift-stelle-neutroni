#ifndef FUN_H
#define FUN_H

/*
 The system to solve is
 dm/dr = f_m = r^3 E(rho)
 dP/dr = f_P = - [m(r) E(rho)/r^2] * [1 + P(r)/E(rho)] * [1 + r^3 P(r)/m(r)] / [1 - 2m(r)/r^2]

 P(rho) = (alpha - 1)a rho^alpha + (beta - 1)b rho^beta
 rho(P) is found by numerically
 */

double f_m(double r, double rho);


double f_P(double r, double rho, double P, double m);


double P_rho(double rho);


double DP_rho(double rho);

double findRho(double P);

#endif
