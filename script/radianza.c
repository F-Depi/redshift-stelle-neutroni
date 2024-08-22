// ::setlocal makeprg=cd\ script\ &&\ gcc\ radianza.c\ -lm\ &&\ ./a.out\ &&\ rm\ a.out
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535     // \pi
#define P0 150.174          // = E_0, MeV/c^2/fm^3
#define R0 20.06145         // km
#define M0 12.655756        // solar masses
#define nu0 2.41798924208e20 //Hz
#define B0 1.30105998011e-6 // Mev / fm^2
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
}


// Metodo Trapezi, usato da <name>_old()
double integrale_trapezio_old(double a, double b, int N, double r, double R, double M, double (*fun)(double, double, double, double)){
    // assume a < b
    double h = (b - a)/ (double)N;
    double integral = ((*fun)(a, r, R, M) + (*fun)(b, r, R, M)) * h / 2;

    for (int i = 1; i < N; i++) {
    integral += h * (*fun)(a + i * h, r, R, M);
    }


    return integral;
}


// Metodo Simpson, usato da <name>_old()
double integrale_simpson_old(double a, double b, int N, double r, double R, double M, double (*fun)(double, double, double, double)){
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
void test_cvg_old(){
    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};        // Raggi delle 3 stelle
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};       // Masse delle 3 stelle
    int N = 10;
    int N_trap, N_simp;
    double Pot, nu_max = 20.;
    double r = R[0];
    double Pot_cvg = integrale_trapezio_old(1e-12, nu_max, 1e8,  r, R[0], M[0], &funB_corrected); // Aggiunto dopo, e' quello ottenuto con N = 1e6
    double errore_max = 1e-5;
    int kk = 0;

    // Dati per grafico cvg per N trapezi
    FILE *f0 = fopen("../data/potenza/test_cvg_N_trap.csv", "w");
    fprintf(f0, "Pot,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_trapezio_old(1e-12, nu_max, N, r, R[0], M[0], &funB_corrected);
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
        Pot = integrale_simpson_old(1e-12, nu_max, N, r, R[0], M[0], &funB_corrected);
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
    Pot_cvg = integrale_trapezio_old(1e-12, 200, 1e8, R[0], R[0], M[0], &funB_corrected);
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

        Pot = integrale_trapezio_old(1e-12, nu_max, N, R[0], R[0], M[0], &funB_corrected);
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

        Pot = integrale_simpson_old(1e-12, nu_max, N, R[0], R[0], M[0], &funB_corrected);
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


// Fissati i parametri per l'integrale con la funzione sopra, calcola il grafico Pot(r) per tutte e 3 le stelle sia con Simpson che con trapezi.
// r_max è in kilometri
void dati_grafico_Pot_old(double N_trap, double N_simp, double nu_max, double r_max){
    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};        // Raggi delle 3 stelle
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};       // Masse delle 3 stelle


    //// Trapezi
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        char filename_trap[50]; sprintf(filename_trap, "../data/potenza/Pot_trap_%d.csv", i + 1);
        FILE *f_trap = fopen(filename_trap, "w");
        fprintf(f_trap, "r,Pot\n");

        double r = R[i];
        double Pot = 0;
        while (r < r_max / R0){
            Pot = integrale_trapezio_old(1e-12, nu_max, N_trap, r, R[i], M[i], &funB_corrected);
            fprintf(f_trap, "%.5e,%.5e\n", r * R0, Pot * POT0);
            r += 0.01;
        }

        // r = \infty
        Pot = integrale_trapezio_old(1e-12, nu_max, N_trap, - 1, R[i], M[i], &funB_corrected);
        fprintf(f_trap, "%.5e,%.5e\n", R[i] * R0, Pot * POT0);

        fclose(f_trap);
    }

    //// Simpson
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        char filename_simp[50]; sprintf(filename_simp, "../data/potenza/Pot_simp_%d.csv", i + 1);
        FILE *f_simp = fopen(filename_simp, "w");
        fprintf(f_simp, "r,Pot\n");

        double r = R[i];
        double Pot = 0;
        while (r < r_max / R0){
            Pot = integrale_simpson_old(1e-12, nu_max, N_simp, r, R[i], M[i], &funB_corrected);
            fprintf(f_simp, "%.5e,%.5e\n", r * R0, Pot * POT0);
            r += 0.01;
        }

        // r = \infty
        Pot = integrale_simpson_old(1e-12, nu_max, N_simp, - 1, R[i], M[i], &funB_corrected);
        fprintf(f_simp, "%.5e,%.5e\n", R[i] * R0, Pot * POT0);

        fclose(f_simp);
    }
}


// Metodo Trapezi
double integrale_trapezio(double a, double b, int N, double T, double (*fun)(double, double)){
    // assume a < b
    double h = (b - a)/ (double)N;
    double integral = ((*fun)(a, T) + (*fun)(b, T)) * h / 2;

    for (int i = 1; i < N; i++) {
    integral += h * (*fun)(a + i * h, T);
    }

    return integral;
}


// Metodo Simpson
double integrale_simpson(double a, double b, int N, double T, double (*fun)(double, double)){
    // assume a < b
    double h = (b - a)/N;
    double integral = ((*fun)(a, T) + (*fun)(b, T)) * h / 3;

    // assume N pari
    if (N % 2 != 0){
        printf("N non è pari");
    }
    for (int i = 1; i <= N/2 - 1; i++) {
        integral += (2 * h /3) * (*fun)(a + 2 * i * h, T);
        integral += (4 * h / 3) * (*fun)(a + (2 * i - 1) * h, T);
    }

    integral += (4 * h / 3 ) * (*fun)(a + (N - 1) * h, T);

    return integral;
}


// Versione semplificata e con solo l'integrale dell'omonima funzione usata sopra
// Viene lasciata la variabile double Pot a indicare il valore dell'integrale solamente
void test_cvg(){
    double T = 0.01;
    int N = 10;
    int N_trap, N_simp;
    double Pot, nu_max = 20.;
    double Pot_cvg = integrale_trapezio(1e-12, nu_max, 1e8, T, &funB);
    double errore_max = 1e-7;
    printf("Errore massimo scelto = %.0e\n", errore_max);
    int kk = 0;

    // Dati per grafico cvg per N trapezi
    FILE *f0 = fopen("../data/potenza/test_cvg_N_trap.csv", "w");
    fprintf(f0, "I,N,nu_max\n");
    while(N <= 1e6){
        Pot = integrale_trapezio(1e-12, nu_max, N, T, &funB);
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
        Pot = integrale_simpson(1e-12, nu_max, N, T, &funB);
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
    T = 1;
    Pot_cvg = integrale_trapezio(1e-12, 200, 1e8, T, &funB);
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

        Pot = integrale_trapezio(1e-12, nu_max, N, T, &funB);
        fprintf(f2, "%.13e,%d,%.13e\n", Pot, N, nu_max);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per trapezi N = %d, nu_max = %.7f sono sufficienti\n", N, nu_max);
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

        Pot = integrale_simpson(1e-12, nu_max, N, T, &funB);
        fprintf(f3, "%.13e,%d,%.13e\n", Pot, N, nu_max);

        if ((fabs(Pot - Pot_cvg) / Pot_cvg) < errore_max && kk == 0){
            printf("Per Simpson N = %d, nu_max = %.7f sono sufficienti\n", N, nu_max);
            kk++;
        }

        nu_max *= 1.05;
    }
    fprintf(f3, "%.13e,%d,%.13e\n", Pot_cvg, (int)1e8, 200.);
    fclose(f3);

}


// Fissati i parametri per l'integrale con la funzione sopra, calcola il grafico Pot(r) per tutte e 3 le stelle sia con Simpson che con trapezi.
// r_max è in kilometri
// Utilizza B senza correzioni, perché è stato fatto un cambio di variabile nell'itegrale
void dati_grafico_Pot(int N_trap, int N_simp, double nu_max, double r_max){

    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};        // Raggi delle 3 stelle
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};       // Masse delle 3 stelle
    double T = 1.;

    double Integrale_trap = integrale_trapezio(1e-12, nu_max, N_trap, T, &funB);
    printf("I con trapezi, N = %d, nu_max = %.1f, vale %.12e\n", N_trap, nu_max, Integrale_trap);
    double corr, Pot;

    //// Trapezi
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        char filename_trap[50]; sprintf(filename_trap, "../data/potenza/Pot_trap_%d.csv", i + 1);
        FILE *f_trap = fopen(filename_trap, "w");
        fprintf(f_trap, "r,Pot\n");

        double r = R[i];
        while(r < r_max / R0){
            corr = pow(1. - 2. * M[i] / R[i], 1. / 2.) * pow(1. - 2. * M[i] / r, - 1. / 2.);
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

    //// Simpson
    double Integrale_simp = integrale_simpson(1e-12, nu_max, N_simp, T, &funB);
    printf("I con Simpson, N = %d, nu_max = %.1f, vale %.12e\n", N_simp, nu_max, Integrale_simp);
    /*
    // Cicla sulle 3 stelle
    for (int i = 0; i < 3; i++){
        char filename_simp[50]; sprintf(filename_simp, "../data/potenza/Pot_simp_%d.csv", i + 1);
        FILE *f_simp = fopen(filename_simp, "w");
        fprintf(f_simp, "r,Pot\n");

        double r = R[i];
        while (r < r_max / R0){
            corr = pow(1. - 2. * M[i] / R[i], 1. / 2.) * pow(1. - 2. * M[i] / r, - 1. / 2.);
            Pot= Integrale_simp * corr;
            fprintf(f_simp, "%.5e,%.5e\n", r * R0, Pot * POT0);
            r += 0.01;
        }

        // r = \infty
        corr = pow(1. - 2. * M[i] / R[i], 1. / 2.);
        Pot = Integrale_simp * corr;
        fprintf(f_simp, "%.5e,%.5e\n", R[i] * R0, Pot * POT0);

        fclose(f_simp);
    }
    */
}


void Teff_su_T(int N_trap, double nu_max, double T_min, double T_max){
    double R[3] = {59.03824 / R0, 10.90280 / R0, 8.559218 / R0};        // Raggi delle 3 stelle
    double M[3] = {14.29963 / M0, 0.9252994 / M0, 1.528782 / M0};       // Masse delle 3 stelle


    // ciclo sulle stelle
    for (int i = 0; i < 3; i++){
        char filename[50]; sprintf(filename, "../data/potenza/Teff_%d.csv", i + 1);
        FILE *f = fopen(filename, "w");
        fprintf(f, "T,Teff\n");

        double T = T_min;
        double Integrale, Teff;

        while (T <= T_max){
            Integrale = integrale_trapezio(1e-12, nu_max, N_trap, T, &funB);
            Teff = pow(15, 1. / 4.) / PI * pow(1. - 2. * M[i] / R[i], 1. / 8.);
            Teff *= pow(Integrale, 1. / 4.);
            fprintf(f, "%.7e,%.7e\n", T, Teff);
            T += 0.001;
        }
        fclose(f);
    }
}


int main(){

    /***************** Radianza corretta *****************/
    // tutto_su_radianza();

    
    /******************* Potenza Totale ******************/
    // test_cvg();

    // int N_trap = 22936;
    // int N_simp = 30516;
    // double nu_max = 24.0661923;
    // double r_max = 500;
    // dati_grafico_Pot(N_trap, N_simp, nu_max, r_max);
    

    /*************** Temperatura Percepita ***************/
    int N_trap = 22936;
    double nu_max = 24.0661923;
    double T_min = 0.01;
    double T_max = 1.;
    Teff_su_T(N_trap, nu_max, T_min, T_max);





    return 0;
}
