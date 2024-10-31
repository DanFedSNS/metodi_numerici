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
const int D = 3;
double beta = 0.46;
int lattice_size;

double min(double x, double y) {
    return (x < y) ? x : y;
}

double dot(double *a, double *b){
    return a[0] * b[0] + a[1] * b[1];
}

double dot_angle(double theta, double phi){
    return cos(theta - phi);
}

double norma(double *a){
    return sqrt(pow(a[0], 2) + pow(a[1], 2));
}

void nearest(int rx, int ry, int rz, int resx[2 * D], int resy[2 * D], int resz[2 * D]) { 
    // Neighbors in the x-direction
    resx[0] = (rx - 1 + L) % L;  // Left neighbor in x
    resx[1] = (rx + 1) % L;      // Right neighbor in x
    resx[2] = rx;                // No change for y-direction neighbor
    resx[3] = rx;                // No change for y-direction neighbor
    resx[4] = rx;                // No change for z-direction neighbor
    resx[5] = rx;                // No change for z-direction neighbor

    // Neighbors in the y-direction
    resy[0] = ry;                // No change for x-direction neighbor
    resy[1] = ry;                // No change for x-direction neighbor
    resy[2] = (ry - 1 + L) % L;  // Down neighbor in y
    resy[3] = (ry + 1) % L;      // Up neighbor in y
    resy[4] = ry;                // No change for z-direction neighbor
    resy[5] = ry;                // No change for z-direction neighbor

    // Neighbors in the z-direction
    resz[0] = rz;                // No change for x-direction neighbor
    resz[1] = rz;                // No change for x-direction neighbor
    resz[2] = rz;                // No change for y-direction neighbor
    resz[3] = rz;                // No change for y-direction neighbor
    resz[4] = (rz - 1 + L) % L;  // Front neighbor in z
    resz[5] = (rz + 1) % L;      // Back neighbor in z
}

int update(double reticolo[L][L][L], int lattice_size){
    bool *reticolo_aus = (bool*) malloc(lattice_size * sizeof(bool));     //reticolo_aus[x][y][z] = reticolo_aus[x * L^2 + y * L + z]
    int *clusterx = (int*) malloc(lattice_size * sizeof(int));
    int *clustery = (int*) malloc(lattice_size * sizeof(int));
    int *clusterz = (int*) malloc(lattice_size * sizeof(int));

    for (int j=0; j<lattice_size; j++){
        reticolo_aus[j] = false;
    }

    int nold = 0;
    int nnew = 1;
    int lc = 1;

    double theta_refl = myrand() * tau;
    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());
    int rz = (int)((double)L * myrand());

    reticolo[rx][ry][rz] = 2 * theta_refl - reticolo[rx][ry][rz];
    clusterx[0] = rx;
    clustery[0] = ry;
    clusterz[0] = rz;
    reticolo_aus[rx * L * L + ry * L + rz] = true;

    int nnx[2*D], nny[2*D], nnz[2*D];
    while (nold < nnew){
        for (int p = nold; p < nnew; p++){
            nearest(clusterx[p], clustery[p], clusterz[p], nnx, nny, nnz);  //nn of cluster[p]
            for (int i = 0; i < 2*D; i++){
                double prob = 1 - exp(min(0.0, 2* beta * (dot_angle(reticolo[clusterx[p]][clustery[p]][clusterz[p]], theta_refl) * dot_angle(reticolo[nnx[i]][nny[i]][nnz[i]], theta_refl))));

                if (reticolo_aus[nnx[i] * L * L + nny[i] * L + nnz[i]] == false && myrand() < prob){
                    reticolo[nnx[i]][nny[i]][nnz[i]] = 2 * theta_refl - reticolo[nnx[i]][nny[i]][nnz[i]];
                    
                    clusterx[lc] = nnx[i];
                    clustery[lc] = nny[i];
                    clusterz[lc] = nnz[i];

                    reticolo_aus[nnx[i] * L * L + nny[i] * L + nnz[i]] = true;
                    lc++;
                }
            }
        }

        nold = nnew;
        nnew = lc;
    }

    free(reticolo_aus);
    free(clusterx);
    free(clustery);
    free(clusterz);

    return lc;
}

double magn(double reticolo[L][L][L]){
    double somma_spin[] = {0.0, 0.0};
    for (int i=0; i < L;  i++){
        for (int j=0; j < L; j++){
            for (int k=0; k < L; k++){
                somma_spin[0] += cos(reticolo[i][j][k]);
                somma_spin[1] += sin(reticolo[i][j][k]);
            }
        }
    }
    
    return norma(somma_spin) / (double) lattice_size;
}

double energy(double reticolo[L][L][L]){    //DA CONTROLLARE
    double sum = 0.0;

    for (int rx = 0; rx < L; rx++) {
        for (int ry = 0; ry < L; ry++) {
            for (int rz = 0; rz < L; rz++) {
                // Interactions in the x-direction
                int rx_next = (rx + 1) % L;
                sum += -dot_angle(reticolo[rx][ry][rz], reticolo[rx_next][ry][rz]);

                // Interactions in the y-direction
                int ry_next = (ry + 1) % L;
                sum += -dot_angle(reticolo[rx][ry][rz], reticolo[rx][ry_next][rz]);

                // Interactions in the z-direction
                int rz_next = (rz + 1) % L;
                sum += -dot_angle(reticolo[rx][ry][rz], reticolo[rx][ry][rz_next]);
            }
        }
    }

    return sum * beta / (double) lattice_size;
}

void initialize_lattice(double lattice[L][L][L]){
    for (int i=0; i < L;  i++){    
        for (int j=0; j < L; j++){
            for (int k=0; k < L; k++){
                lattice[i][j][k] = myrand() * tau;
            }
        } 
    }
}

void termalizzazione(double lattice[L][L][L], int iterations, int lattice_size){
    printf("\nAggiornamenti sulla termalizzazione:\n");
    int lunghezza_cluster_media = 0;
    for (int i=0; i < iterations; i++){  //termalizzazione
        lunghezza_cluster_media += update(lattice, lattice_size); 
        
        if (i % 1000 == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)i, (double)iterations);
        }
    }

    printf("\nlunghezza media cluster = %f", (double)lunghezza_cluster_media / iterations);
}

void presa_misure(double lattice[L][L][L], int lattice_size, int num_measures, int iter_bet_meas, FILE *fp, FILE *fp_config, bool save_config){
    printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas; i++){
            update(lattice, lattice_size);
        }

        fprintf(fp, "%f, %f\n", magn(lattice), energy(lattice));    //questo rallenta molto
        
        if (save_config == true){
            for (int a=0; a < L; a++){
                for (int b=0; b < L; b++){
                    for (int c=0; c < L; c++){
                        fprintf(fp_config, "%f, ", lattice[a][b][c]);
                    }
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

    double lattice[L][L][L];    //potrebbe essere meglio usare malloc (e restrict)
    lattice_size = pow(L, D);

    initialize_lattice(lattice);
    
    char datafile[40], datafile_config[50]; // file name
    FILE *fp, *fp_config; // pointer to file

    sprintf(datafile, "./xy_model_cluster/L%d_beta%.2f.dat", L, beta); // file name initialized with a string
    fp = fopen(datafile, "w");
    fprintf(fp, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);    
    
    if (save_config == true){
        sprintf(datafile_config, "./config/xy_model_cluster_L%d_beta%.2f.dat", L, beta);
        fp_config = fopen(datafile_config, "w");
        fprintf(fp_config, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);
    }

    termalizzazione(lattice, iterations, lattice_size);
    
    presa_misure(lattice, lattice_size, num_measures, iter_bet_meas, fp, fp_config, save_config);


    fclose(fp);
    if (save_config == true){
        fclose(fp_config);    
    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nHa impiegato %.3f secondi\n", time_spent);


    char datafile_time[] = "./time/xy_model_cluster.dat";
    FILE *fp_time;

    fp_time = fopen(datafile_time, "a");
    fprintf(fp_time, "\ntime = %.2f min, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", time_spent/60, beta, L, iterations, iter_bet_meas, num_measures);    
    fclose(fp_time);
}


int main(void){
    int iterations = 1e6;
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