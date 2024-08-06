// ::setlocal makeprg=cd\ script\ &&\ gcc\ punto7.c\ -lm\ &&\ ./a.out\ &&\ rm\ a.out
#include <stdio.h>
#include <math.h>
#define PI 3.1415926535     // \pi
#define P0 150.174          // = E_0, MeV/c^2/fm^3
#define R0 20.06145         // km
#define M0 12.655756        // solar masses
#define nu0 2.41798924208e20 //Hz
#define B0 1.30105998011e-6 // Mev / fm^2
#define POT0 (PI * B0 * nu0)


// Legge i dati generati nei 3 file data/MR_N.csv, serve solo la lunghezza del file -1 come input
int read_MR_data(int tipo_politropica, int lenfile, double *P, double *R, double *M){

    char fdata_name[50];
    sprintf(fdata_name, "../data/MR_%d.csv", tipo_politropica);

    // Open file where the data computed above is stored
    FILE *file = fopen(fdata_name, "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    
    // Read P, R, M
    double temp;
    if (fscanf(file, "%c,%c0,%c,%c", &temp, &temp, &temp, &temp) != 4){
            printf("Error reading header from file\n");
            return 1;
    }
    for (int i = 0; i < lenfile; i++) {
        double col1, col2, col3;
        if (fscanf(file, "%*lf,%lf,%lf,%lf", &col1, &col2, &col3) != 3){
            printf("i = %d\n", i);
            printf("Error reading data from file\n");
            return 1;
        }
        P[i] = col1 / P0;
        R[i] = col2 / R0;
        M[i] = col3 / M0;
    }
    
    fclose(file);

    return 0;
}

// Radianza
double funB(double nu, double T){
    return pow(nu, 3) / (exp(nu / T) - 1.);
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


double funTeff(double R, double M, double integral_T){
    return pow(15., 1. / 4.) / PI * pow(1. - 2. * M / R, 1. / 8.) * pow(integral_T, 1. / 4.);
}

int main(){

    int file_lens[3] = {395, 617, 717};
    int N_trap = 271;
    double nu_max = 30.;
    double T;

    double integral_T1 = integrale_trapezio(1e-12, nu_max, N_trap, 1.,  &funB);

    for (int i = 0; i < 3; i++){

        double P[file_lens[i]], R[file_lens[i]], M[file_lens[i]];
        read_MR_data(i, file_lens[i], P, R, M);

        char filename[50]; sprintf(filename, "../data/punto7/Teff1_su_P_%d.csv", i + 1);
        FILE *f = fopen(filename, "w");
        fprintf(f, "P,Teff\n");

        for (int j = 0; j < file_lens[i]; j++){
            T = funTeff(R[j], M[j], integral_T1);
            fprintf(f, "%.7e,%.7e\n", P[j] * P0, T);
        }
        fclose(f);
    }

    double integral_T001 = integrale_trapezio(1e-12, nu_max, N_trap*10, 0.01,  &funB);

    for (int i = 0; i < 3; i++){

        double P[file_lens[i]], R[file_lens[i]], M[file_lens[i]];
        read_MR_data(i, file_lens[i], P, R, M);

        char filename[50]; sprintf(filename, "../data/punto7/Teff001_su_P_%d.csv", i + 1);
        FILE *f = fopen(filename, "w");
        fprintf(f, "P,Teff\n");

        for (int j = 0; j < file_lens[i]; j++){
            T = funTeff(R[j], M[j], integral_T001);
            fprintf(f, "%.7e,%.7e\n", P[j] * P0, T);
        }
        fclose(f);
    }

    return 0;
}
