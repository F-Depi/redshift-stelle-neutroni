#include <stdio.h>
#include <math.h>
#include "fun.h"
#define P_0 150.174       // = E_0, MeV/fm^3
#define R_0 20.06145      // km
#define M_0 12.655756     // solar masses
#define a (13.4 / P_0)      // energy density parameter
#define b (5.6 / P_0)       // energy density parameter
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
    printf("a = %f\n", a);
    printf("b = %f\n", b);
    printf("alpha = %f\n", alpha);
    printf("beta = %f\n", beta);
    printf("alpha - 1 = %f\n", alpha1);
    printf("beta - 1 = %f\n", beta1);
    printf("alpha1 * a = %f\n", alpha1 * a);
    printf("beta1 * b = %f\n", beta1 * b);


    double rho, P;
    FILE *f1 = fopen("../data/rho-P.csv", "w");
    fprintf(f1, "rho,P\n");
    for (int i = 0; i < 1000; i++){
        rho = i*0.01;
        fprintf(f1, "%e,%e\n", rho, P_rho(rho));
    }

    FILE *f2 = fopen("../data/P-rho.csv", "w");
    fprintf(f2, "P,rho\n");
    for (int i = 0; i < 25000; i++){
        P = i*0.01;
        fprintf(f2, "%e,%e\n", P, findRho(P));
    }
}













