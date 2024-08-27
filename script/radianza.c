// ::setlocal makeprg=cd\ script\ &&\ gcc\ radianza.c\ -lm\ &&\ ./a.out\ &&\ rm\ a.out
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535             // \pi
#define P0 150.33046048             // = E_0, MeV/fm^3
#define R0 19.996542543             // km
#define M0 13.542058427             // solar masses
#define nu0 2.41798924208e20        //Hz
#define B0 1.30105998011e-6         // Mev / fm^2
#define POT0 (PI * B0 * nu0)


// Radianza
double funB(double nu, double T){
    return pow(nu, 3) / (exp(nu / T) - 1.);
}


// Radianza corretta con redshift gravitazionale
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


// Esercizio sulla Radianza
void tutto_su_radianza(){

    // Calcolo una semplice radianza
    double T = 1.;
    int N = 1000;
    double B, nu, nu_max = 15.;

    FILE *f1 = fopen("../data/radianza/B.csv", "w");
    fprintf(f1, "nu,B\n");
    for (int i = 0; i < N; i++){
        nu = nu_max * (double)i / (double)N;
        B = funB(nu, T);
        fprintf(f1, "%.10e,%.10e\n", nu * nu0, funB(nu, T) * B0);
    }
    fclose(f1);


    // Prendo i dati di Phi dal file (a mano)
    double R[3] = {11.04289 / R0, 10.86752 / R0, 8.531525 / R0};
    double M[3] = {2.456841 / M0, 0.990100 / M0, 1.635845 / M0};
    double r[3] = {1.5, 8., -1};
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        // Cicla sulle 3 distanze
        for (int j = 0; j < 3; j++){

            char filename3[50];
            sprintf(filename3, "../data/radianza/Bcorr_%d_%.1f.csv", i + 1, r[j]);
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
}


// Metodo Trapezi
double integrale_trapezio(double a, double b, int N, double T,
                                                double (*fun)(double, double)){
    // assume a < b
    double h = (b - a)/ (double)N;
    double integral = ((*fun)(a, T) + (*fun)(b, T)) * h / 2;

    for (int i = 1; i < N; i++) {
    integral += h * (*fun)(a + i * h, T);
    }

    return integral;
}


// Metodo Simpson
double integrale_simpson(double a, double b, int N, double T,
                                                double (*fun)(double, double)){
    // assume a < b
    double h = (b - a)/N;
    double integral = ((*fun)(a, T) + (*fun)(b, T)) * h / 3;

    // assume N pari
    if (N % 2 != 0){
        printf("N non e' pari");
    }
    for (int i = 1; i <= N/2 - 1; i++) {
        integral += (2 * h /3) * (*fun)(a + 2 * i * h, T);
        integral += (4 * h / 3) * (*fun)(a + (2 * i - 1) * h, T);
    }

    integral += (4 * h / 3 ) * (*fun)(a + (N - 1) * h, T);

    return integral;
}


// Verifica della convergenza dell'integrale della radianza per i metodi simpson
// e trapezio. Prima viene fissato l'estremo superiore e fatto variare N,
// trovato l'N sufficiente per l'errore voluto si fissa N/nu_max e si fa variare
// nu_max fino a raggiungere l'errore voluto.
// Si prende come riferimento l'integrale calcolato con N = 1e8 nu_max = 20
// per la prima parte. N = 1e8 nu_max = 200 per la seconda parte.
void test_cvg(){
    int N = 10;
    int N_trap, N_simp;
    double Pot, nu_max = 20.;
    double Pot_cvg = integrale_trapezio(1e-12, nu_max, 1e8, 1., &funB);
    double errore_max = 1e-7;
    printf("Errore massimo scelto = %.0e\n", errore_max);
    int kk = 0;

    // Dati per grafico cvg per N trapezi
    FILE *f0 = fopen("../data/potenza/test_cvg_N_trap.csv", "w");
    fprintf(f0, "I,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_trapezio(1e-12, nu_max, N, 1., &funB);
        fprintf(f0, "%.13e,%d,%.13e\n", Pot, N, nu_max * nu0);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per trapezi N = %d e' sufficiente\n", N);
            N_trap = N;
            kk++;
        }

        N *= 1.1;
        if (N % 2 != 0) N += 1;
    }
    fprintf(f0, "%.13e,%d,%.13e\n", Pot_cvg, (int)1e8, nu_max * nu0);
    fclose(f0);


    N = 10;
    kk = 0;
    // Dati per grafico cvg per N simpson
    FILE *f1 = fopen("../data/potenza/test_cvg_N_simp.csv", "w");
    fprintf(f1, "Pot,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_simpson(1e-12, nu_max, N, 1., &funB);
        fprintf(f1, "%.13e,%d,%.13e\n", Pot, N, nu_max * nu0);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per Simpson N = %d e' sufficiente\n", N);
            N_simp = N;
            kk++;
        }

        N *= 1.1;
        if (N % 2 != 0) N += 1;
    }
    fprintf(f1, "%.13e,%d,%.13e\n", Pot_cvg, (int)1e8, nu_max * nu0);
    fclose(f1);

    // Decidiamo che va bene N_trap e N_simp risultati del codice sopra
    // Ora verifichiamo la convergenza di nu_max nel caso peggiore, ovvero
    Pot_cvg = integrale_trapezio(1e-12, 200, 1e8, 1., &funB);
    nu_max = 20.;
    int Nrel_trap = (double)N_trap / nu_max; // Teniamo la stessa densita'
    int Nrel_simp = (double)N_simp / nu_max; // di trapezi: N / nu_max

    // partiamo da
    nu_max = 10.;
    kk = 0;
    // Test cvg per A (ovvero nu_max)
    FILE *f2 = fopen("../data/potenza/test_cvg_A_trap.csv", "w");
    fprintf(f2, "Pot,N,nu_max[ad]\n");
    while(nu_max <= 1e2){
        N = Nrel_trap * nu_max;
        if (N % 2 != 0) N += 1;

        Pot = integrale_trapezio(1e-12, nu_max, N, 1., &funB);
        fprintf(f2, "%.13e,%d,%.13e\n", Pot, N, nu_max);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per trapezi N = %d, nu_max = %.7f sono sufficienti\n",
                                                                    N, nu_max);
            kk++;
        }

        nu_max *= 1.05;
    }
    fprintf(f2, "%.13e,%d,%.13e\n", Pot_cvg, (int)1e8, 200.);
    fclose(f2);

    nu_max = 10.;
    kk = 0;
    // Test cvg per A (ovvero nu_max)
    FILE *f3 = fopen("../data/potenza/test_cvg_A_simp.csv", "w");
    fprintf(f2, "Pot,N,nu_max[ad]\n");
    while(nu_max <= 1e2){
        N = Nrel_simp * nu_max;
        if (N % 2 != 0) N += 1;

        Pot = integrale_simpson(1e-12, nu_max, N, 1., &funB);
        fprintf(f3, "%.13e,%d,%.13e\n", Pot, N, nu_max);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per Simpson N = %d, nu_max = %.7f sono sufficienti\n",
                                                                    N, nu_max);
            kk++;
        }

        nu_max *= 1.05;
    }
    fprintf(f3, "%.13e,%d,%.13e\n", Pot_cvg, (int)1e8, 200.);
    fclose(f3);
}


// Fissati i parametri per l'integrale con la funzione sopra, calcola il grafico
// Pot(r) per tutte e 3 le stelle con i trapezi.
// r_max è in kilometri
// Utilizza B senza correzioni, perché è stato fatto un cambio di variabile
// nell'itegrale
void dati_grafico_Pot(int N_trap, double nu_max, double r_max){

    double R[3] = {11.04289 / R0, 10.86752 / R0, 8.531525 / R0}; // Raggi stelle
    double M[3] = {2.456841 / M0, 0.990100 / M0, 1.635845 / M0}; // Masse stelle
    double T = 1.;

    double Integrale_trap = pow(T, 4) * integrale_trapezio(1e-12, nu_max,
                                                            N_trap, 1., &funB);

    printf("I con trapezi, N = %d, nu_max = %.1f, vale %.12e\n",
            N_trap, nu_max, Integrale_trap);
    double corr, Pot;

    //// Trapezi
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        char filename_trap[50];
        sprintf(filename_trap, "../data/potenza/Pot_trap_%d.csv", i + 1);
        FILE *f_trap = fopen(filename_trap, "w");
        fprintf(f_trap, "r,Pot\n");

        double r = R[i];
        while(r < r_max / R0){
            corr = pow(1. - 2. * M[i] / R[i], 1. / 2.)
                 * pow(1. - 2. * M[i] / r, - 1. / 2.);
            Pot = Integrale_trap * corr;
            fprintf(f_trap, "%.5e,%.5e\n", r * R0, Pot * POT0);
            r += 0.01;
        }

        // r = \infty
        corr = pow(1. - 2. * M[i] / R[i], 1. / 2.);
        Pot = Integrale_trap * corr;
        fprintf(f_trap, "%.5e,%.5e\n", R[i] * R0, Pot * POT0);

        fclose(f_trap);
    }
}


void Teff_su_T(int N_trap, double nu_max, double T_min, double T_max){
    double R[3] = {11.04289 / R0, 10.86752 / R0, 8.531525 / R0}; // Raggi stelle
    double M[3] = {2.456841 / M0, 0.990100 / M0, 1.635845 / M0}; // Masse stelle


    // ciclo sulle stelle
    for (int i = 0; i < 3; i++){
        char filename[50];
        sprintf(filename, "../data/potenza/Teff_%d.csv", i + 1);
        FILE *f = fopen(filename, "w");
        fprintf(f, "T,Teff\n");

        double T = T_min;
        double Teff, I = integrale_trapezio(1e-12, nu_max, N_trap, 1., &funB);
        double I4 = pow(I, 1. / 4.);

        while (T <= T_max){
            Teff = pow(15, 1. / 4.) / PI * pow(1. - 2. * M[i] / R[i], 1. / 8.);
            Teff *= T * I4;
            fprintf(f, "%.7e,%.7e\n", T, Teff);
            T += 0.001;
        }
        fclose(f);
    }
}


int main(){

    /***************** Radianza corretta *****************/

    //tutto_su_radianza();

    
    /******************* Potenza Totale ******************/
    
    //test_cvg();

    //int N_trap = 264;
    //double nu_max = 29.2526072;
    //double r_max = 500;
    //dati_grafico_Pot(N_trap, nu_max, r_max);
    

    /*************** Temperatura Efficacie ***************/
    
    //int N_trap = 264;
    //double nu_max = 29.2526072;
    //double T_min = 0.001;
    //double T_max = 1.;
    //Teff_su_T(N_trap, nu_max, T_min, T_max);


    return 0;
}
