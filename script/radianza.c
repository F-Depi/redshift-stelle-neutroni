// :setlocal makeprg=cd\ script\ &&\ gcc\ radianza.c\ -lm\ &&\ ./a.out\ &&\ rm\ a.out
#include <stdio.h>
#include <math.h>
#define P0 150.174          // = E_0, MeV/c^2/fm^3
#define R0 20.06145         // km
#define M0 12.655756        // solar masses
#define nu0 2.41798924208e20 //Hz
#define A 0.00161310878747 // Mev / fm^2


double funB(double nu){
    return pow(nu, 3) / (exp(nu) - 1.);
}

double funB_corrected(double nu, double r, double R, double M){
    double corr;
    if (r < 0)
        corr = pow(1. - 2 * M / R, - 1. / 2);
    else
        corr = pow((1. - 2 * M / r) / (1. - 2 * M / R), 1. / 2);
    nu *= corr;
    return pow(nu, 3) / (exp(nu) - 1.);
}


int main(){

    // Calcolo una semplice radianza
    int N = 1000;
    double B, nu, nu_max = 15.;

    FILE *f1 = fopen("../data/radianza/B.csv", "w");
    fprintf(f1, "nu,B\n");
    for (int i = 0; i < N; i++){
        nu = nu_max * (double)i / (double)N;
        B = funB(nu);
        fprintf(f1, "%.10e,%.10e\n", nu * nu0, funB(nu) * A);
    }
    fclose(f1);


    /****************** Radianza corretta ********************/
    // Prendo i dati di Phi dal file

    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};
    double r[3] = {1.5, 8., -1};

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            char filename[50]; sprintf(filename, "../data/radianza/Bcorr_%d_%.1f.csv", i + 1, r[j]);
            FILE *f = fopen(filename, "w");
            fprintf(f, "nu,B\n");
            for (int k = 0; k < N; k++){
                nu = nu_max * (double)k / (double)N;
                B = funB_corrected(nu, r[j] * R[i], R[i], M[i]);
                fprintf(f, "%.10e,%.10e\n", nu * nu0, B * A);
            }
            fclose(f);
        }
    }
    return 0;
}
