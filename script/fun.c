#include <stdio.h>
#include <math.h>
#include "fun.h"
#define P0 150.174          // = E_0, MeV/c^2/fm^3
#define R0 20.06145         // km
#define M0 12.655756        // solar masses
#define RHO0 0.16           // fm^-3
#define A (13.4 / P0)       // energy density parameters
#define B (5.6 / P0)       
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


// f_m = r^3 E
double fun_m(double r, double P, int tipo_politropica){
    return r * r * fun_E(P, tipo_politropica);
}


// f_P = - (P + E)(m + r^3 P)/(r^2 - 2mr)
double fun_P(double r, double P, double m, int tipo_politropica){
    if (m == 0)
        return 0;
    return (P + fun_E(P, tipo_politropica)) * (m + pow(r, 3) * P) / ((2 * m - r) * r);
}


// This is for findRho()
double P_of_rho(double rho){
    return ALPHA1 * A * pow(rho, ALPHA) + BETA1 * B * pow(rho, BETA);
}


// This is for findRho(), it's the analitic derivative of the function above
double DP_of_rho(double rho){
    return ALPHA * ALPHA1 * A * pow(rho, ALPHA - 1) + BETA * BETA1 * B * pow(rho, BETA - 1);
}


// Since E = E(rho), we need to find rho from P to evaluate E
double findRho(double P){

    // Newton-Raphson method to find rho
    if (P > 0.01){

        double rho = pow(P / (BETA1 * B), 1 / BETA);    // good approximation
        while (fabs(P - P_of_rho(rho)) > 1e-6){
            rho -= (P_of_rho(rho) - P) / DP_of_rho(rho);
        }
        return rho;
    }

    // for small P it's safer to use the bisect method
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


// Evaluates the energy for 3 differents matter state functions
// tipo_politropica = 1, 2, 3
double fun_E(double P, int tipo_politropica){
    // Politropica quasi realistica (a*rho^alpha + b*rho^beta)
    if (tipo_politropica == 0){
        double rho = findRho(P);
        return A * pow(rho, ALPHA) + B * pow(rho, BETA);
    }

    double lambda, K;

    // Materia fermionica non relativistica
    if (tipo_politropica == 1){
        lambda = 5. / 3.;
        K = 0.05;
    } else if (tipo_politropica == 2){
        lambda = 2.54;
        K = 0.01;
    } else {
        printf("Tipo politropica non riconosciuto\n");
        return 0;
    }

    double a1 = P / (lambda - 1);
    return a1 + pow(a1 / K, 1. / lambda);
}


//RK4 algorithm to advance 1 step of length h, in a sistem like dP/dr = fun_P dm/dr = fun_m 
// Updates the values of m, P given
// tipo_politropica = 1, 2, 3
void rungeKutta4(double h, double r, double *P, double *m, int tipo_politropica){

    double k1, k2, k3, k4, l1, l2, l3, l4;

    k1 = h * fun_m(r, *P, tipo_politropica);
    l1 = h * fun_P(r, *P, *m, tipo_politropica);

    k2 = h * fun_m(r + h / 2, *P + l1 / 2, tipo_politropica);
    l2 = h * fun_P(r + h / 2, *P + l1 / 2, *m + k1 / 2, tipo_politropica);

    k3 = h * fun_m(r + h / 2, *P + l2 / 2, tipo_politropica);
    l3 = h * fun_P(r + h / 2, *P + l2 / 2, *m + k2 / 2, tipo_politropica);

    k4 = h * fun_m(r + h, *P + l3, tipo_politropica);
    l4 = h * fun_P(r + h, *P + l3, *m + k3, tipo_politropica);

    *m += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    *P += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
}



