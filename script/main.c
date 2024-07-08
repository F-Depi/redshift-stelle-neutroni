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


void solve_system(double h, double r, double P, double m, double RM[2], int tipo_politropica){
    while (P > 0){
        RM[0] = r;
        RM[1] = m;
        r += h;
        rungeKutta4(h, r, &P, &m, tipo_politropica);
        if (r < 2 * m){
            printf("BH");
            break;
        }
    }
}


void get_MR(int tipo_politropica){

    char filename[50];
    sprintf(filename, "../data/MR_%d.csv", tipo_politropica);
    FILE *f = fopen(filename, "w");
    fprintf(f, "h,P0,R,M\n");

    double startP, P, prevM, m, r;
    double RM[2] = {};
    double h = 1e-5;

    int start_i;
    if (tipo_politropica == 0){
        start_i = -450;
    } else if (tipo_politropica == 1){
        start_i = -600;
    } else if (tipo_politropica == 2){
        start_i = -650;
    }

    for (int i = start_i; i < 500; i++){

        prevM = RM[1];
        startP = pow(10, i / 100.);
        P = startP;
        m = 0;
        r = 0;
        RM[0] = 0;
        RM[1] = 0;

        solve_system(h, r, P, m, RM, tipo_politropica);

        if (RM[1] < prevM){
            printf("Limit reached at P = %.2e MeV c^-3 fm^-3\n", startP * P0);
            break;
        }
        if (RM[0] * R0 < 3){
            printf("Star too small, not saving the data\n");
            continue;
        }

        fprintf(f, "%e,%e,%e,%e\n", h, startP * P0, RM[0] * R0, RM[1] * M0);
        fflush(f);
        printf("i = %d, P = %.2e\n", i, startP);
    }

    fclose(f);

}


int main(){

    #pragma omp parallel for num_threads(3)
    for (int i = 1; i < 3; i++){
        get_MR(i);
    }

    return 0;
}













