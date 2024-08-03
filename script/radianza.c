// ::setlocal makeprg=cd\ script\ &&\ gcc\ radianza.c\ -lm\ &&\ ./a.out\ &&\ rm\ a.out
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535     // \pi
#define P0 150.174          // = E_0, MeV/c^2/fm^3
#define R0 20.06145         // km
#define M0 12.655756        // solar masses
#define nu0 2.41798924208e20 //Hz
#define B0 1.30105998011e-6 // Mev / fm^2
#define POT0 (2 * PI * B0 * nu0)


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
    for (int i = 1; i <= N/2 - 1; i++) {
        integral += (2 * h /3) * (*fun)(a + 2 * i * h, r, R, M);
        integral += (4 * h / 3) * (*fun)(a + (2 * i - 1) * h, r, R, M);
    }


    integral += (4 * h / 3 ) * (*fun)(a + (N - 1) * h, r, R, M);

    return integral;
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

            char filename[50]; sprintf(filename, "../data/radianza/Bcorr_%d_%.1f.csv", i + 1, r[j]);
            FILE *f = fopen(filename, "w");
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
    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};
    double r[3] = {1.5, 8., -1};
    int N = 10;
    double Pot, nu_max = 20.;

    // Test cvg per N trapezi
    FILE *f0 = fopen("../data/potenza/test_cvg_N_trap.csv", "w");
    fprintf(f0, "Pot,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_trapezio(1e-10, nu_max, N, r[0] * R[0], R[0], M[0], &funB_corrected);
        fprintf(f0, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max * nu0);
        N *= 1.1;
    }
    fclose(f0);

    N = 10;

    // Test cvg per N simpson
    FILE *f1 = fopen("../data/potenza/test_cvg_N_simp.csv", "w");
    fprintf(f1, "Pot,N,nu_max\n");
    while(2 * N <= 1e6){
        Pot = integrale_simpson(1e-10, nu_max, 2 * N, r[0] * R[0], R[0], M[0], &funB_corrected);
        fprintf(f1, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max * nu0);
        N *= 1.1;
    }
    fclose(f1);

    // Decidiamo che va bene
    nu_max = 10;
    N = 1e3;
    int Nrel = (double)N / nu_max;
    // Test cvg per A (ovvero nu_max)
    FILE *f2 = fopen("../data/potenza/test_cvg_A_trap.csv", "w");
    fprintf(f2, "Pot,N,nu_max[ad]\n");
    while(nu_max <= 1e2){
        N = Nrel * nu_max;
        Pot = integrale_trapezio(1e-10, nu_max, N, r[0] * R[0], R[0], M[0], &funB_corrected);
        fprintf(f2, "%.13e,%d,%.13e\n", Pot * POT0, N, nu_max);
        nu_max *= 1.1;
    }
    fclose(f2);

/*
// Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        // Cicla sulle 3 distanze
        for (int j = 0; j < 3; j++){

            char filename[50]; sprintf(filename, "../data/potenza/Pot_trap_%d_%.1f.csv", i + 1, r[j]);
            FILE *f = fopen(filename, "w");
            fprintf(f, "nu,Pot\n");

            for (int k = 0; k < N; k++){
                nu = nu_max * (double)k / (double)N;
                B = funB_corrected(nu, r[j] * R[i], R[i], M[i]);
                fprintf(f, "%.10e,%.10e\n", nu * nu0, B * A);
            }
            fclose(f);
        }
    }
    */

    return 0;
}
