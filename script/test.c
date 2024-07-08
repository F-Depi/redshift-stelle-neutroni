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

 P(rho) = (alpha - 1)a rho^alpha + (beta - 1)b rho^beta
 rho(P) is found by numerically
 */


void test_cnts_findRho(){
    printf("a = %f\n", A);
    printf("b = %f\n", B);
    printf("alpha = %f\n", ALPHA);
    printf("beta = %f\n", BETA);
    printf("alpha - 1 = %f\n", ALPHA1);
    printf("beta - 1 = %f\n", BETA1);
    printf("alpha1 * a = %f\n", ALPHA1 * A);
    printf("beta1 * b = %f\n", BETA1 * B);


    double rho, P;
    FILE *f1 = fopen("../data/rho-P.csv", "w");
    fprintf(f1, "rho,P\n");
    for (int i = 0; i < 1000; i++){
        rho = i*0.01;
        fprintf(f1, "%e,%e\n", rho * R0, P_of_rho(rho) * P0);
    }

    FILE *f2 = fopen("../data/P-rho.csv", "w");
    fprintf(f2, "P,rho\n");
    for (int i = 0; i < 25000; i++){
        P = i*0.01;
        fprintf(f2, "%e,%e\n", P * P0, findRho(P) * RHO0);
    }

}


void test_cvg(){
    double m = 0;
    double P = 1;
    double r = 0;
    double rho = findRho(P);
    double h = 1e-5;


    FILE *f = fopen("../data/data2.csv", "w");
    fprintf(f, "r,P,m,rho\n");
    fprintf(f, "%e,%e,%e,%e\n", r * R0, P * P0, m * M0, rho * RHO0);

    while (P > 0){           // Analizzando la curva di P si vede che non si scende di molto sotto 0.01 (si arriva circa a 0.004), quindi 0.01 va bene.
        r += h;
        rungeKutta4(h, r, &P, &m, 2);
        rho = findRho(P);
        fprintf(f, "%e,%e,%e,%e\n", r * R0, P * P0, m * M0, rho * RHO0);
        fflush(f);
    }

    fclose(f);
}


void test_h(){
    FILE *f2 = fopen("../data/data_cvg_1.csv", "w");
    fprintf(f2, "h,R,M,P\n");

    for (double h = 1e-2; h > 1e-8; h /= 2){
        printf("h = %e\n", h);
        double R_stella = -3.14;
        double M_stella = -3.14;
        double P_supercifie = -3.14;
        double m = 0;
        double P = P0;
        double r = 0;

        while (P > 0){
            // Saving data before the last step, to avoid the case where P < 0
            R_stella = r;
            M_stella = m;
            P_supercifie = P;

            r += h;
            rungeKutta4(h, r, &P, &m, 1);
        }
        fprintf(f2, "%e,%e,%e,%e\n", h, R_stella * R0, M_stella * M0, P_supercifie * P0);
    }

    fclose(f2);
}


int main(){

    // Test findRho
    //test_cnts_findRho();

    
    // Funziona? Converge?
    //test_cvg();
    // Sì


    // What's the best h increment?
    //test_h();
    // h = 1e-4 should be fine, but let's use 1e-5 to be sure
}













