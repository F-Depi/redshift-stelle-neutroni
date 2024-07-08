#include <stdio.h>
#include <math.h>
#include <omp.h>
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


void solve_system(double h, double r, double P, double m, double RM[2]){
    while (P > 0){
        RM[0] = r;
        RM[1] = m;
        r += h;
        rungeKutta4(h, r, &P, &m);
        if (r < 2 * m){
            printf("BH");
            break;
        }
    }
}


int main(){

    FILE *f1 = fopen("../data/RM_medium.csv", "w");
    fprintf(f1, "h,P0,R,M\n");

    #pragma omp parallel for schedule(dynamic), shared(f1), num_threads(8)
    for (int i = -100; i < 100; i++){
        double startP = pow(10, i / 10.);
        double P = startP;
        double RM[2] = {};
        double m = 0;
        double r = 0;
        double h = 1e-5;

        solve_system(h, r, P, m, RM);

        #pragma omp critical
        {
            fprintf(f1, "%e,%e,%e,%e\n", h, startP * P0, RM[0] * R0, RM[1] * M0);
            fflush(f1);
            printf("i = %d, P = %.2e\n", i, startP);
        }
    }

    fclose(f1);

    return 0;
}













