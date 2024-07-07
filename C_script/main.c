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



int main(){

    double m = 0.00001;
    double P = P0;
    double r = 0;
    double rho = findRho(P);
    double h = 0.0001;
    FILE *f = fopen("../data/data.csv", "w");
    fprintf(f, "r,P,m,rho\n");
    fprintf(f, "%e,%e,%e,%e\n", r, P, m, rho);

    while (P > 0.1){
        r += h;
        rungeKutta4(h, r, &P, &m);
        rho = findRho(P);
        fprintf(f, "%e,%e,%e,%e\n", r, P, m, rho);
    }

    fclose(f);

    return 0;
}













