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
    fprintf(f, "h,P0,R,lenfile\n");

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


// Salva in 3 file diversi r,P,m di 3 stelle con 3 politropiche diverse, date le 3 pressioni iniziali
void gen_maxM_data(double *startP){

    #pragma omp parallel for num_threads(3)
    for (int i = 0; i < 3; i++){
    int tipo_politropica = i;
    double h = 1e-5;
    double r = 0;
    double m = 0;
    double P = startP[i];

    char fdata_name[50];
    sprintf(fdata_name, "../data/maxM_%d.csv", tipo_politropica);
    FILE *f = fopen(fdata_name, "w");
    fprintf(f, "r,P,m\n");

    while (P > 0){
        fprintf(f, "%.10e,%.10e,%.10e\n", r * R0, P * P0, m * M0);
        r += h;
        rungeKutta4(h, r, &P, &m, tipo_politropica);
    }
    fclose(f);

    }
}



// Legge i dati generati nei 3 file sopra, serve solo la lunghezza del file -1 come input
int read_maxM_data(int tipo_politropica, int lenfile, double *r, double *P, double *m){

    char fdata_name[50];
    sprintf(fdata_name, "../data/maxM_%d.csv", tipo_politropica);

    // Open file where the data computed above is stored
    FILE *file = fopen(fdata_name, "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    
    // Read R, m
    double temp;
    if (fscanf(file, "%c,%c,%c", &temp, &temp, &temp) != 3) {
            printf("Error reading header from file\n");
            return 1;
    }
    for (int i = 0; i < lenfile; i++) {
        double col1, col2, col3;
        if (fscanf(file, "%lf,%lf,%lf", &col1, &col2, &col3) != 3) {
            printf("i = %d\n", i);
            printf("Error reading data from file\n");
            return 1;
        }
        r[i] = col1 / R0;
        P[i] = col2 / P0;
        m[i] = col3 / M0;
    }
    
    fclose(file); // close file

    return 0;
}

int main(){

    /* Genera i dati per costruire il grafico massa raggio, richiesto dal punto 2 */
    // #pragma omp parallel for num_threads(3)
    // for (int i = 0; i < 3; i++){
    //     get_MR(i);
    // }


    /* Genera dati per fare Grafico del potenziale, richiesto dal punto 3 */
    // double startP[3] = {43.31065 / P0, 217.0675 / P0, 947.5339 / P0};
    // gen_maxM_data(startP);



    /* Legge i dati dai 3 file generati sopra */
    int len_files[3] = {294288, 54348, 42666};
    int tipo_politropica = 0;

    int lenfile = len_files[tipo_politropica];
    double r[lenfile] = {};
    double P[lenfile] = {};
    double m[lenfile] = {};

    read_maxM_data(tipo_politropica, lenfile, r, P, m);

    return 0;
}













