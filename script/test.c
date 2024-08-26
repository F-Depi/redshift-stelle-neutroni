// ::setlocal makeprg=cd\ script\ &&\ make\ test\ &&\ ./test.x
#include <stdio.h>
#include "fun.h"

/*
 The system to solve is
 dm/dr = f_m = r^2 E(rho)
 dP/dr = f_P = - (P + E)(m + r^3 P)/(r^2 - 2mr)

 P(rho) = ALPHA A rho^(ALPHA+1) + BETA B rho^(BETA+1)
 rho(P) is found numerically
*/


void print_constants(){
    printf("a = %f\n", A);
    printf("b = %f\n", B);
    printf("alpha = %f\n", ALPHA);
    printf("alpha + 1 = %f\n", ALPHA + 1);
    printf("beta = %f\n", BETA);
    printf("alpha * a = %f\n", ALPHA * A);
    printf("beta * b = %f\n", BETA * B);
    printf("b beta /a /alhpa = %f\n", 5.62 * BETA / 13.4 / ALPHA);
    printf("P(rho = 0.16 fm^-3) = %f\n", P_of_rho(1.));
}


void test_rhovsP(){
    double rho = 0.001;
    double P = 0;
    FILE *f1 = fopen("../data/test/rho_of_P.csv", "w");
    fprintf(f1, "rho,P\n");
    while (P < 1){
        P = P_of_rho(rho);
        fprintf(f1, "%e,%e\n", rho * N0, P_of_rho(rho) * P0);
        rho += 0.001;
    }

    P = 0.0001;
    FILE *f2 = fopen("../data/test/P_of_rho.csv", "w");
    fprintf(f2, "P,rho\n");
    while(P < 1){
        fprintf(f2, "%e,%e\n", P * P0, findRho(P) * N0);
        P *= 1.1;
    }

}


void compare_eneries(){
    for (int i = 0; i < 3; i++){
        char filename[50]; sprintf(filename, "../data/test/E_of_P_%d.csv", i + 1);
        FILE *f = fopen(filename, "w");
        fprintf(f, "rho,P\n");

        double P = 0;
        while(P <= 1){
            fprintf(f, "%e,%e\n", P * P0, fun_E(P, i) * P0);
            P += 0.001;
        }
        fclose(f);
    }
}


void test_cvg_stella(double startP, double smallestP, int tipo_politropica){
    double m = 0;
    double P = startP;
    double r = 0;
    double rho = findRho(P);
    double h = 1e-5;


    char filename[50]; sprintf(filename, "../data/test/data%d.csv", tipo_politropica + 1);
    FILE *f = fopen(filename, "w");
    fprintf(f, "r,P,m,rho\n");
    fprintf(f, "%e,%e,%e,%e\n", r * R0, P * P0, m * M0, rho * N0);

    while (P > smallestP){
        fprintf(f, "%e,%e,%e,%e\n", r * R0, P * P0, m * M0, rho * N0);
        r += h;
        rungeKutta4(h, r, &P, &m, tipo_politropica);
        double prevP = P;
        rho = findRho(P);
        if (r < 2 * m){
            printf("BH");
            break;
        }
        if (P > prevP){
            printf("sei arrivato in superficie\n");
            break;
        }
    }

    fclose(f);
}


void test_cvg_h(){
    FILE *f = fopen("../data/test/cvg_1.csv", "w");
    fprintf(f, "h,R,M,P\n");

    for (double h = 1e-2; h > 1e-8; h /= 2){
        printf("h = %e\n", h);
        double R_stella = -3.14;
        double M_stella = -3.14;
        double P_supercifie = -3.14;
        double m = 0;
        double P = 0.1;
        double r = 0;

        while (P > 0){
            // Saving data before the last step, to avoid the case where P < 0
            R_stella = r;
            M_stella = m;
            P_supercifie = P;

            r += h;
            rungeKutta4(h, r, &P, &m, 0);
        }
        fprintf(f, "%e,%e,%e,%e\n", h, R_stella * R0, M_stella * M0, P_supercifie * P0);
    }

    fclose(f);
}


int main(){

    /************** Test findRho **************/
    // Test is the energy equations plus related are correct
    // Test if the constant are written properly
    // Test if the results at least appears self consistent
    // print_constants();
    // test_rhovsP();
    compare_eneries();

    
    /************** First simulation *************/
    // Tries to solve the equations once for 1 star, to see if
    // everything is ok and the pressure converges to 0.
    // Take the politropic type as input
    // test_cvg_stella(2, 0.00001, 0);
    // test_cvg_stella(2., 0., 1);
    // test_cvg_stella(2., 0., 2);


    /*************** Test h ***********************/
    // Tries different increments to see how that affects the results
    // only for 1 star, it shouldn't matter
    // test_cvg_h();
    // h = 1e-4 should be fine, but let's use 1e-5 to be sure
}













