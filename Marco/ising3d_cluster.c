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
#define pos(rx, ry, rz, L) (rx * L * L + ry * L + rz)


int L = 10;
const int D = 3;
double beta = 0.23;
int lattice_size;
char modello[] = "ising3d_cluster";

int update(int *restrict reticolo, int lattice_size, double prob){
    bool *restrict reticolo_aus = (bool*) malloc(lattice_size * sizeof(bool));     //reticolo_aus[x][y][z] = reticolo_aus[x * L^2 + y * L + z]
    int *restrict clusterx = (int*) malloc(lattice_size * sizeof(int));
    int *restrict clustery = (int*) malloc(lattice_size * sizeof(int));
    int *restrict clusterz = (int*) malloc(lattice_size * sizeof(int));

    for (int j=0; j<lattice_size; j++){
        reticolo_aus[j] = false;
    }

    int nold = 0;
    int nnew = 1;
    int lc = 1;

    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());
    int rz = (int)((double)L * myrand());

    clusterx[0] = rx;
    clustery[0] = ry;
    clusterz[0] = rz;
    reticolo_aus[pos(rx, ry, rz, L)] = true;

    int nnx[2*D], nny[2*D], nnz[2*D];
    while (nold < nnew){
        for (int p = nold; p < nnew; p++){
            nearest_cu(clusterx[p], clustery[p], clusterz[p], nnx, nny, nnz, L);  //nn of cluster[p]
            for (int i = 0; i < 2*D; i++){
                if (reticolo_aus[pos(nnx[i], nny[i], nnz[i], L)] == false && reticolo[pos(nnx[i], nny[i], nnz[i], L)] == reticolo[pos(clusterx[p], clustery[p], clusterz[p], L)] && myrand() < prob){                    
                    clusterx[lc] = nnx[i];
                    clustery[lc] = nny[i];
                    clusterz[lc] = nnz[i];

                    reticolo_aus[pos(nnx[i], nny[i], nnz[i], L)] = true;
                    lc++;
                }
            }
        }

        nold = nnew;
        nnew = lc;
    }

    for (int i = 0; i < lc; i++){
        reticolo[pos(clusterx[i], clustery[i], clusterz[i], L)] *= -1;
    }

    free(reticolo_aus);
    free(clusterx);
    free(clustery);
    free(clusterz);

    return lc;
}

double magn(int *restrict reticolo, int lattice_size){
    int somma_spin = 0;
    for (int i=0; i < lattice_size;  i++){
            somma_spin += reticolo[i];
    }
    
    return (double) somma_spin / (double) lattice_size;
}

double energy(int *restrict reticolo, int lattice_size){    //DA CONTROLLARE
    double sum = 0.0;

    for (int rx = 0; rx < L; rx++) {
        for (int ry = 0; ry < L; ry++) {
            for (int rz = 0; rz < L; rz++) {
                // Interactions in the x-direction
                int rx_next = (rx + 1) % L;
                sum += -reticolo[pos(rx, ry, rz, L)] * reticolo[pos(rx_next, ry, rz, L)];

                // Interactions in the y-direction
                int ry_next = (ry + 1) % L;
                sum += -reticolo[pos(rx, ry, rz, L)] * reticolo[pos(rx, ry_next, rz, L)];

                // Interactions in the z-direction
                int rz_next = (rz + 1) % L;
                sum += -reticolo[pos(rx, ry, rz, L)] * reticolo[pos(rx, ry, rz_next, L)];
            }
        }
    }

    return sum * beta / (double) lattice_size;
}

void initialize_lattice(int *restrict lattice, int lattice_size){
    for (int i=0; i < lattice_size;  i++){    
        lattice[i] = 1;
    }
}  

void termalizzazione(int *restrict lattice, int iterations, int lattice_size, double prob){
    printf("\nAggiornamenti sulla termalizzazione:\n");
    int lunghezza_cluster_media = 0;
    for (int i=0; i < iterations; i++){  //termalizzazione
        lunghezza_cluster_media += update(lattice, lattice_size, prob); 
        
        if (i % 1000 == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)i, (double)iterations);
        }
    }

    printf("\nlunghezza media cluster = %f", (double)lunghezza_cluster_media / iterations);
}

void presa_misure(int *restrict lattice, int lattice_size, int num_measures, int iter_bet_meas, FILE *restrict fp, FILE *restrict fp_config, bool save_config, double prob){
    printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas; i++){
            update(lattice, lattice_size, prob);
        }

        fprintf(fp, "%f, %f\n", magn(lattice, lattice_size), energy(lattice, lattice_size));    //questo rallenta molto
        
        if (save_config == true){
            for (int a=0; a < lattice_size; a++){
                fprintf(fp_config, "%d, ", lattice[a]);
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

    lattice_size = pow(L, D);
    int *restrict lattice = (int *) malloc(lattice_size * sizeof(int));
    double prob = 1 -  exp(- 2 * beta);

    initialize_lattice(lattice, lattice_size);
    
    FILE *fp, *fp_config; //pointer to file

    init_file(modello, L, beta, &fp, iterations, iter_bet_meas, num_measures, save_config, &fp_config);

    termalizzazione(lattice, iterations, lattice_size, prob);
    
    presa_misure(lattice, lattice_size, num_measures, iter_bet_meas, fp, fp_config, save_config, prob);

    free(lattice);

    close_file(&fp, save_config, &fp_config);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nHa impiegato %.3f secondi\n", time_spent);

    save_time_spent(beta, L, modello, time_spent, iterations, iter_bet_meas, num_measures);
}


int main(void){
    int iterations = 1e5;
    int iter_bet_meas = 10;    //iterations between two measures
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