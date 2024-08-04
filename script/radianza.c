// ::setlocal makeprg=cd\ script\ &&\ gcc\ radianza.c\ -lm\ &&\ ./a.out\ &&\ rm\ a.out
#include <complex.h>
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535     // \pi
#define P0 150.174          // = E_0, MeV/c^2/fm^3
#define R0 20.06145         // km
#define M0 12.655756        // solar masses
#define nu0 2.41798924208e20 //Hz
#define B0 1.30105998011e-6 // Mev / fm^2
#define POT0 (PI * B0 * nu0)


double funB(double nu){
    return pow(nu, 3) / (exp(nu) - 1.);
}

double funB_corrected(double nu, double r, double R, double M){
    double corr;
    // To simulate r = infinity
    if (r < 0)
        corr = pow(1. - 2 * M / R, - 1. / 2);
    else
        corr = pow((1. - 2 * M / r) / (1. - 2 * M / R), 1. / 2);
    nu *= corr;
    return pow(nu, 3) / (exp(nu) - 1.);
}


// Metodo Trapezi
double integrale_trapezio(double a, double b, int N, double r, double R, double M, double (*fun)(double, double, double, double)){
    // assume a < b
    double h = (b - a)/ (double)N;
    double integral = ((*fun)(a, r, R, M) + (*fun)(b, r, R, M)) * h / 2;

    for (int i = 1; i < N; i++) {
    integral += h * (*fun)(a + i * h, r, R, M);
    }


    return integral;
}


// Simpson method to integrate
double integrale_simpson(double a, double b, int N, double r, double R, double M, double (*fun)(double, double, double, double)){
    // assume a > b
    double h = (b - a)/N;
    double integral = ((*fun)(a, r, R, M) + (*fun)(b, r, R, M)) * h / 3;

    // assume N pari
    if (N % 2 != 0){
        printf("N non è pari");
    }
    for (int i = 1; i <= N/2 - 1; i++) {
        integral += (2 * h /3) * (*fun)(a + 2 * i * h, r, R, M);
        integral += (4 * h / 3) * (*fun)(a + (2 * i - 1) * h, r, R, M);
    }


    integral += (4 * h / 3 ) * (*fun)(a + (N - 1) * h, r, R, M);

    return integral;
}


// Test di convergenza ripetto ai parametri nu_max e N
void test_cvg(){
    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};        // Raggi delle 3 stelle
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};       // Masse delle 3 stelle
    int N = 10;
    int N_trap, N_simp;
    double Pot, nu_max = 20.;
    double r = R[0];
    double Pot_cvg = integrale_trapezio(1e-10, nu_max, 1e8,  r, R[0], M[0], &funB_corrected); // Aggiunto dopo, e' quello ottenuto con N = 1e6
    double errore_max = 1e-5;
    int kk = 0;

    // Dati per grafico cvg per N trapezi
    FILE *f0 = fopen("../data/potenza/test_cvg_N_trap.csv", "w");
    fprintf(f0, "Pot,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_trapezio(1e-10, nu_max, N, r, R[0], M[0], &funB_corrected);
        fprintf(f0, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max * nu0);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per trapezi N = %d e' sufficiente\n", N);
            N_trap = N;
            kk++;
        }

        N *= 1.1;
        if (N % 2 != 0) N += 1;
    }
    fprintf(f0, "%.13e,%d,%.13e\n", Pot_cvg * POT0, (int)1e8, nu_max * nu0);
    fclose(f0);


    N = 10;
    kk = 0;
    // Dati per grafico cvg per N simpson
    FILE *f1 = fopen("../data/potenza/test_cvg_N_simp.csv", "w");
    fprintf(f1, "Pot,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_simpson(1e-12, nu_max, N, r, R[0], M[0], &funB_corrected);
        fprintf(f1, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max * nu0);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per Simpson N = %d e' sufficiente\n", N);
            N_simp = N;
            kk++;
        }

        N *= 1.1;
        if (N % 2 != 0) N += 1;
    }
    fprintf(f1, "%.13e,%d,%.13e\n", Pot_cvg * POT0, (int)1e8, nu_max * nu0);
    fclose(f1);

    // Decidiamo che va bene N_trap e N_simp risultati del codice sopra
    Pot_cvg = integrale_trapezio(1e-12, 200, 1e8, R[0], R[0], M[0], &funB_corrected);
    nu_max = 20.;
    int Nrel_trap = (double)N_trap / nu_max; // Teniamo la stella densita'
    int Nrel_simp = (double)N_simp / nu_max; // di N / nu_max

    // partiamo da
    nu_max = 10.;
    kk = 0;
    // Test cvg per A (ovvero nu_max)
    FILE *f2 = fopen("../data/potenza/test_cvg_A_trap.csv", "w");
    fprintf(f2, "Pot,N,nu_max[ad]\n");
    while(nu_max <= 1e2){
        N = Nrel_trap * nu_max;
        if (N % 2 != 0) N += 1;

        Pot = integrale_trapezio(1e-10, nu_max, N, R[0], R[0], M[0], &funB_corrected);
        fprintf(f2, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per trapezi N = %d, nu_max = %.7f sono sufficienti\n", N, nu_max);
            kk++;
        }

        nu_max *= 1.05;
    }
    fprintf(f2, "%.13e,%d,%.13e\n", Pot_cvg * POT0, (int)1e8, 200.);
    fclose(f2);

    nu_max = 10.;
    kk = 0;
    // Test cvg per A (ovvero nu_max)
    FILE *f3 = fopen("../data/potenza/test_cvg_A_simp.csv", "w");
    fprintf(f2, "Pot,N,nu_max[ad]\n");
    while(nu_max <= 1e2){
        N = Nrel_simp * nu_max;
        if (N % 2 != 0) N += 1;

        Pot = integrale_simpson(1e-10, nu_max, N, R[0], R[0], M[0], &funB_corrected);
        fprintf(f3, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per Simpson N = %d, nu_max = %.7f sono sufficienti\n", N, nu_max);
            kk++;
        }

        nu_max *= 1.05;
    }
    fprintf(f3, "%.13e,%d,%.13e\n", Pot_cvg * POT0, (int)1e8, 200.);
    fclose(f3);

}

int main(){

    /*
    // Calcolo una semplice radianza
    int N = 1000;
    double B, nu, nu_max = 15.;

    FILE *f1 = fopen("../data/radianza/B.csv", "w");
    fprintf(f1, "nu,B\n");
    for (int i = 0; i < N; i++){
        nu = nu_max * (double)i / (double)N;
        B = funB(nu);
        fprintf(f1, "%.10e,%.10e\n", nu * nu0, funB(nu) * B0);
    }
    fclose(f1);
    */


    /****************** Radianza corretta ********************/
    /*
    // Prendo i dati di Phi dal file (a mano)
    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};
    double r[3] = {1.5, 8., -1};
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        // Cicla sulle 3 distanze
        for (int j = 0; j < 3; j++){

            char filename3[50]; sprintf(filename3, "../data/radianza/Bcorr_%d_%.1f.csv", i + 1, r[j]);
            FILE *f = fopen(filename3, "w");
            fprintf(f, "nu,B\n");
            for (int k = 0; k < N; k++){
                nu = nu_max * (double)k / (double)N;
                B = funB_corrected(nu, r[j] * R[i], R[i], M[i]);
                fprintf(f, "%.10e,%.10e\n", nu * nu0, B * B0);
            }
            fclose(f);
        }
    }
    */

    
    /****************** Potenza Totale ********************/
    test_cvg();

    // Actual exercise
    //N = 1e3;
    //nu_max = 20.;

    ////// Trapezi
    //// Cicla sulle 3 stelle
    //char filename3[50]; sprintf(filename3, "../data/potenza/Pot_simp_%d.csv", i);
    //FILE *f3 = fopen(filename3, "w");
    //fprintf(f3, "r,R,M,Pot,Int,N,nu_max[ad]\n");

    //for (int i = 0; i < 3; i++){
    //    // Cicla sulle 3 distanze
    //    for (int j = 0; j < 3; j++){

    //        Pot = integrale_trapezio(1e-10, nu_max, N, r[j] * R[i], R[i], M[i], &funB_corrected);
    //        fprintf(f3, "%.10e,%.10e,%.10e,%.10e,%s,%d,%.0f\n", r[j] * R[i] * R0, R[i] * R0, M[i] * M0, Pot * POT0, "simpson", N, nu_max);
    //        Pot = integrale_simpson(1e-10, nu_max, N, r[j] * R[i], R[i], M[i], &funB_corrected);
    //        fprintf(f3, "%.10e,%.10e,%.10e,%.10e,%s,%d,%.0f\n", r[j] * R[i] * R0, R[i] * R0, M[i] * M0, Pot * POT0, "trapezi", N / 2, nu_max);
    //    }
    //}
    //fclose(f3);


    return 0;
}
