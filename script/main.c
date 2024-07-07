#include <stdio.h>
#include <math.h>
#include <omp.h>
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

    FILE *f1 = fopen("../data/RM2.csv", "w");
    fprintf(f1, "h,P0,R,M\n");

    #pragma omp parallel for schedule(dynamic), shared(f1), num_threads(8)
    for (int i = 0; i < 80; i++){
        double P = (pow(10 , i/ 8.) - 0.9) * P0;
        printf("i = %d, P = %.2e\n", i, P);
        double R_stella = 0;
        double M_stella = 0;
        double m = 0;
        double r = 0;
        double h = 1e-5;

        while (P > 0.01){
            r += h;
            rungeKutta4(h, r, &P, &m);
            R_stella = r;
            M_stella = m;
        }
        #pragma omp critical
        {
            fprintf(f1, "%e,%e,%e,%e\n", h, P0 / (double)i, R_stella, M_stella);
            fflush(f1);
        }
    }

    fclose(f1);

    return 0;
}













