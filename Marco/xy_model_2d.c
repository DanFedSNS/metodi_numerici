#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "./include/random.h"
#include "./include/get_array.h"
#include <omp.h>

#define pi 3.14159265358979323846
#define tau 6.28318530717958647692
//each spin is an angle theta with components (cos theta, sin theta)

int L = 10;
const int D = 2;
double beta = .47;
int lattice_size;
double metro_or_micro = 0.15;


double dot(double *a, double *b){
    return a[0] * b[0] + a[1] * b[1];
}

double dot_angle(double theta, double phi){
    return cos(theta - phi);
}

double norma(double *a){
    return sqrt(pow(a[0], 2) + pow(a[1], 2));
}

void nearest(int rx, int ry, int resx[2*D], int resy[2*D]){  //restituisce i nearest neighbours di r
    resx[0] = (rx - 1 + L) % L;     //+ L serve a evitare negativi
    resx[1] = (rx + 1) % L;
    resx[2] = rx;
    resx[3] = rx;

    resy[0] = ry;
    resy[1] = ry;
    resy[2] = (ry - 1 + L) % L;
    resy[3] = (ry + 1) % L;
}

void update_metro(double reticolo[L][L]){
    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());

    int nnx[2*D], nny[2*D];
    nearest(rx, ry, nnx, nny);

    double new_spin = myrand() * tau;
    double deltaE = 0; //beta * (E' - E)
    for (int i = 0; i < 2*D; i++){
        deltaE += -dot_angle(new_spin, reticolo[nnx[i]][nny[i]]) + dot_angle(reticolo[rx][ry], reticolo[nnx[i]][nny[i]]);
    }
    deltaE *= beta;

    if (deltaE <= 0){
        reticolo[rx][ry] = new_spin;
    }

    else {
        double w = myrand();
        if (w <= exp(-deltaE)){
            reticolo[rx][ry] = new_spin;
        }
    }
}

void update_micro(double reticolo[L][L]){
    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());

    int nnx[2*D], nny[2*D];
    nearest(rx, ry, nnx, nny);

    double Sr[] = {0.0, 0.0};
    for (int i = 0; i < 2*D; i++){
        Sr[0] += cos(reticolo[nnx[i]][nny[i]]);
        Sr[1] += sin(reticolo[nnx[i]][nny[i]]);
    }

    double Sr_angle = atan2(Sr[1], Sr[0]);  //angolo tra Sr e x-axis, compreso tra -pi e pi, non è un problema

    reticolo[rx][ry] = 2 * Sr_angle - reticolo[rx][ry]; // Riflessione, theta -> 2 * Sr_angle - theta
}

void update(double reticolo[L][L]){
    double r_mm = myrand(); //random per decidere metro or micro
        
    if (r_mm < metro_or_micro){
        for (int j = 0; j < lattice_size; j++){
            update_metro(reticolo);
            }
    }

    else{
        for (int j = 0; j < lattice_size; j++){
            update_micro(reticolo);
            }
    }
}

double magn(double reticolo[L][L]){
    double somma_spin[] = {0.0, 0.0};
    for (int i=0; i < L;  i++){
        for (int j=0; j < L; j++){
            somma_spin[0] += cos(reticolo[i][j]);
            somma_spin[1] += sin(reticolo[i][j]);
        }
    }
    
    return norma(somma_spin) / (double) lattice_size;
}

double energy(double reticolo[L][L]){    //DA CONTROLLARE
    double sum = 0.0;

    for (int rx = 0; rx < L; rx++) {
        for (int ry = 0; ry < L; ry++) {
            int rx_next = (rx + 1) % L;
            sum += -dot_angle(reticolo[rx][ry], reticolo[rx_next][ry]);

            int ry_next = (ry + 1) % L;
            sum += -dot_angle(reticolo[rx][ry], reticolo[rx][ry_next]);
        }
    }

    return sum * beta / (double) lattice_size;
}

void initialize_lattice(double lattice[L][L]){
    for (int i=0; i < L;  i++){   
        for (int j=0; j < L; j++){
            lattice[i][j] = myrand() * tau;
        } 
    }
}

void termalizzazione(double lattice[L][L], int iterations){
    printf("\nAggiornamenti sulla termalizzazione:\n");
    for (int i=0; i < iterations; i++){  //termalizzazione
        update(lattice); 
        
        if (i % 1000 == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)i, (double)iterations);
        }
    }
}

void presa_misure(double lattice[L][L], int num_measures, int iter_bet_meas, FILE *fp, bool save_config, FILE *fp_config){
    printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas; i++){
            update(lattice);
        }

        fprintf(fp, "%f, %f\n", magn(lattice), energy(lattice));
        
        if (save_config == true){
            for (int a=0; a < L; a++){
                for (int b=0; b < L; b++){
                    fprintf(fp_config, "%f, ", lattice[a][b]);
                }
            }       
            fseek(fp_config, -2, SEEK_CUR); //rimuove ", " finali
            fprintf(fp_config, "\n");
        }

        if (j % 1000 == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)j, (double)num_measures);
        }
    }

}

void montecarlo(int iterations, int iter_bet_meas, int num_measures, bool save_config){
    clock_t begin = clock();

    double lattice[L][L];    //potrebbe essere meglio usare malloc (e restrict)
    lattice_size = pow(L, D);

    initialize_lattice(lattice);
    
    char datafile[40], datafile_config[40]; // file name
    FILE *fp, *fp_config; // pointer to file

    sprintf(datafile, "./xy_model_2d/L%d_beta%.2f.dat", L, beta); // file name initialized with a string
    fp = fopen(datafile, "w");
    fprintf(fp, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);    
    
    if (save_config == true){
        sprintf(datafile_config, "./config_xy/L%d_beta%.2f.dat", L, beta);
        fp_config = fopen(datafile_config, "w");
        fprintf(fp_config, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);
    }
    
    termalizzazione(lattice, iterations);    
    
    presa_misure(lattice, num_measures, iter_bet_meas, fp, save_config, fp_config);

    fclose(fp);
    if (save_config == true){
        fclose(fp_config);    
    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nHa impiegato %.3f secondi\n", time_spent);


    char datafile_time[] = "./time/xy_model_2d.dat";
    FILE *fp_time;

    fp_time = fopen(datafile_time, "a");
    fprintf(fp_time, "\ntime = %.2f min, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", time_spent/60, beta, L, iterations, iter_bet_meas, num_measures);    
    fclose(fp_time);
}


int main(void){
    int iterations = 1e5;
    int iter_bet_meas = 1;    //iterations between two measures
    int num_measures = 1e4;
    bool save_config = false;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;
    myrand_init(seed1, seed2);

    /*
    //scegli i valori di beta e L (ed eventualmente degli altri parametri)
    const char *datafile_Larray, *datafile_betaarray;
    datafile_Larray = "L_array.txt";
    datafile_betaarray = "beta_array.txt";

    int num_L = count_lines(datafile_Larray);
    int num_beta = count_lines(datafile_betaarray);

    int L_array[num_L];
    double beta_array[num_beta];

    get_array_from_txt_int(datafile_Larray, L_array);
    get_array_from_txt_double(datafile_betaarray, beta_array);

    //omp_set_num_threads(2);
    //#pragma omp parallel for collapse(2)  bisogna passare  beta e L come input
    for (int i = 0; i < num_beta; i++){
        beta = beta_array[i];
        for (int j = 0; j < num_L; j++){
            L = L_array[j];

            metropolis();
        }
    }
    */

    montecarlo(iterations, iter_bet_meas, num_measures, save_config);

    return EXIT_SUCCESS;
}