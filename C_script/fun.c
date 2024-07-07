#include <stdio.h>
#include <math.h>
#include "fun.h"
#define P0 150.174       // = E_0, MeV/fm^3
#define R0 20.06145      // km
#define M0 12.655756     // solar masses
#define a (13.4 / P0)      // energy density parameter
#define b (5.6 / P0)       // energy density parameter
#define alpha 0.514
#define beta 3.436
#define alpha1 (alpha - 1)
#define beta1 (beta - 1)

/*
 The system to solve is
 dm/dr = f_m = r^3 E(rho)
 dP/dr = f_P = - [m(r) E(rho)/r^2] * [1 + P(r)/E(rho)] * [1 + r^3 P(r)/m(r)] / [1 - 2m(r)/r^2]

 P(rho) = (alpha - 1)a rho^alpha + (beta - 1)b rho^beta
 rho(P) is found by numerically
 */

double f_m(double r, double P){
    double rho = findRho(P);
    double rho3 = pow(rho, 3);
    return pow(r, 3) * ( a * rho3 + b * rho);
}


double f_P(double r, double P, double m){
    if (m == 0)
        return 0;
    double rho = findRho(P);
    double rho3 = pow(rho, 3);
    double Er = a * rho3 + b * rho;
    return - m * Er * (1 + P / Er) * (1 + pow(r, 3) * P / m) / (r - 2 * m);
}


double P_rho(double rho){
    return alpha1 * a * pow(rho, alpha) + beta1 * b * pow(rho, beta);
}


double DP_rho(double rho){
    return alpha * alpha1 * a * pow(rho, alpha - 1) + beta * beta1 * b * pow(rho, beta - 1);
}

double findRho(double P){
    // Newton-Raphson method to find rho

    double rho = pow(P / (beta1 * b), 1 / beta);    // good approximation

    while (fabs(P - P_rho(rho)) > 1e-6){
        rho -= (P_rho(rho) - P) / DP_rho(rho);
    }
    return rho;
}

void rungeKutta4(double h, double r, double *P, double *m){

    double k1, k2, k3, k4, l1, l2, l3, l4;

    k1 = h * f_m(r, *P);
    l1 = h * f_P(r, *P, *m);

    k2 = h * f_m(r + h / 2, *P + l1 / 2);
    l2 = h * f_P(r + h / 2, *P + l1 / 2, *m + k1 / 2);

    k3 = h * f_m(r + h / 2, *P + l2 / 2);
    l3 = h * f_P(r + h / 2, *P + l2 / 2, *m + k2 / 2);

    k4 = h * f_m(r + h, *P + l3);
    l4 = h * f_P(r + h, *P + l3, *m + k3);

    *m += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    *P += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
}
