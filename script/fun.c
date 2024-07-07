#include <stdio.h>
#include <math.h>
#include "fun.h"
#define P0 150.174       // = E_0, MeV/fm^3
#define R0 20.06145      // km
#define M0 12.655756     // solar masses
#define A (13.4 / P0)      // energy density parameter
#define B (5.6 / P0)       // energy density parameter
#define ALPHA 0.514
#define BETA 3.436
#define ALPHA1 (ALPHA - 1)
#define BETA1 (BETA - 1)

/*
 The system to solve is
 dm/dr = f_m = r^3 E(rho)
 dP/dr = f_P = - [m(r) E(rho)/r^2] * [1 + P(r)/E(rho)] * [1 + r^3 P(r)/m(r)] / [1 - 2m(r)/r^2]

 P(rho) = (ALPHA - 1)A rho^ALPHA + (BETA - 1)B rho^BETA
 rho(P) is found by numerically
 */

double fun_m(double r, double P){
    return pow(r, 3) * fun_E(P);
}


double fun_P(double r, double P, double m){
    if (m == 0)
        return 0;
    return (fun_E(P) + P) * (m + pow(r, 3) * P) / ((2 * m - r) * r);
}


double P_of_rho(double rho){
    return ALPHA1 * A * pow(rho, ALPHA) + BETA1 * B * pow(rho, BETA);
}


double DP_of_rho(double rho){
    return ALPHA * ALPHA1 * A * pow(rho, ALPHA - 1) + BETA * BETA1 * B * pow(rho, BETA - 1);
}

double findRho(double P){

    // Newton-Raphson method to find rho
    if (P > 0.01){

        double rho = pow(P / (BETA1 * B), 1 / BETA);    // good approximation
        while (fabs(P - P_of_rho(rho)) > 1e-6){
            rho -= (P_of_rho(rho) - P) / DP_of_rho(rho);
        }
        return rho;
    }

    // Bisect method
    double rho_sx = 0.7, rho_dx = 1.;
    double rho = (rho_sx + rho_dx) / 2;
    while (fabs(P - P_of_rho(rho)) > 1e-6){
        if (P_of_rho(rho) > P)
            rho_dx = rho;
        else
            rho_sx = rho;
        rho = (rho_sx + rho_dx) / 2;
    }
    return rho;
}

double fun_E(double P){
    double rho = findRho(P);
    return A * pow(rho, ALPHA) + B * pow(rho, BETA);
}


void rungeKutta4(double h, double r, double *P, double *m){

    double k1, k2, k3, k4, l1, l2, l3, l4;

    k1 = h * fun_m(r, *P);
    l1 = h * fun_P(r, *P, *m);

    k2 = h * fun_m(r + h / 2, *P + l1 / 2);
    l2 = h * fun_P(r + h / 2, *P + l1 / 2, *m + k1 / 2);

    k3 = h * fun_m(r + h / 2, *P + l2 / 2);
    l3 = h * fun_P(r + h / 2, *P + l2 / 2, *m + k2 / 2);

    k4 = h * fun_m(r + h, *P + l3);
    l4 = h * fun_P(r + h, *P + l3, *m + k3);

    *m += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    *P += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
}
