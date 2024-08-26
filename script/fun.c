#include <stdio.h>
#include <math.h>
#include "fun.h"

/*
 The system to solve is
 dm/dr = f_m = r^2 E(rho)
 dP/dr = f_P = - (P + E)(m + r^3 P)/(r^2 - 2mr)

 P(rho) = ALPHA A rho^(ALPHA+1) + BETA B rho^(BETA+1)
 rho(P) is found numerically
*/


// f_m = r^2 E
double fun_m(double r, double P, int tipo_politropica){
    return r * r * fun_E(P, tipo_politropica);
}


// f_P = - (P + E)(m + r^3 P)/(r^2 - 2mr)
double fun_P(double r, double P, double m, int tipo_politropica){
    if (m == 0)
        return 0;
    return (P + fun_E(P, tipo_politropica)) * (m + pow(r, 3) * P)
           / ((2 * m - r) * r);
}


// This is used in findRho() to do Newton-Rapson
double P_of_rho(double rho){
    return ALPHA * A * pow(rho, ALPHA + 1.) + BETA * B * pow(rho, BETA + 1.);
}


// This is used in findRho() to do Newton-Rapson
// It's the analitic derivative of the function above
double DP_of_rho(double rho){
    return ALPHA * (ALPHA + 1.) * A * pow(rho, ALPHA)
           + BETA * (BETA + 1.) * B * pow(rho, BETA);
}


// Since E = E(rho), we need to invert numerically P(rho) to find rho given a P
double findRho(double P){

    // Newton-Raphson method to find rho
    double rho = pow(P / 0.15, 1. / 3.);    // decent approximation for P
                                           
    while (fabs(P - P_of_rho(rho)) > 1e-8){
        rho -= (P_of_rho(rho) - P) / DP_of_rho(rho);
    }
    return rho;
}


// Evaluates the energy for 3 differents matter state functions
// tipo_politropica = 1, 2, 3
double fun_E(double P, int tipo_politropica){

    // Politropica quasi realistica (rho*mc^2 + a*rho^alpha + b*rho^beta)
    if (tipo_politropica == 0){
        double rho = findRho(P);
        return rho + A * pow(rho, ALPHA + 1.) + B * pow(rho, BETA + 1.);
    }

    double lambda, K;

    // Materia fermionica non relativistica
    if (tipo_politropica == 1){
        lambda = 5. / 3.;
        K = 0.05;
    }
    else if (tipo_politropica == 2){
        lambda = 2.54;
        K = 0.01;
    }
    else {
        printf("Tipo politropica non riconosciuto\n");
        return 0;
    }

    double a1 = P / (lambda - 1.);
    return a1 + pow(a1 / K, 1. / lambda);
}


// RK4 algorithm to advance 1 step of length h, in a sistem like
// dP/dr = fun_P
// dm/dr = fun_m 
// Updates the values of m, P given
// tipo_politropica = 1, 2, 3
void rungeKutta4(double h, double r, double *P, double *m, int tipo_politropica){

    double k1, k2, k3, k4, l1, l2, l3, l4;

    k1 = h * fun_m(r, *P, tipo_politropica);
    l1 = h * fun_P(r, *P, *m, tipo_politropica);

    k2 = h * fun_m(r + h / 2., *P + l1 / 2., tipo_politropica);
    l2 = h * fun_P(r + h / 2., *P + l1 / 2., *m + k1 / 2., tipo_politropica);

    k3 = h * fun_m(r + h / 2., *P + l2 / 2., tipo_politropica);
    l3 = h * fun_P(r + h / 2., *P + l2 / 2., *m + k2 / 2., tipo_politropica);

    k4 = h * fun_m(r + h, *P + l3, tipo_politropica);
    l4 = h * fun_P(r + h, *P + l3, *m + k3, tipo_politropica);

    *m += (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
    *P += (l1 + 2. * l2 + 2. * l3 + l4) / 6.;
}



