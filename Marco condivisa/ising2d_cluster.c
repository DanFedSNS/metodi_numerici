#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "./include/random.h"
#include "./include/get_array.h"
#include <omp.h>
#define pos(rx, ry, L) (rx * L + ry)

int L;
const int D = 2;
int q;
double beta;
int lattice_size;
void (*nearest)(int, int, int *, int *, int);
double (*energy)(int *restrict, int, int, double);
char modello[] = "ising2d_hex_cluster"; //ising2d_tri_cluster or ising2d_sq_cluster or ising2d_hex_cluster 

int update(int *restrict reticolo, int lattice_size, double prob){
    bool *restrict reticolo_aus = (bool*) malloc(lattice_size * sizeof(bool));     //reticolo_aus[x][y][z] = reticolo_aus[x * L^2 + y * L + z]
    int *restrict clusterx = (int*) malloc(lattice_size * sizeof(int));
    int *restrict clustery = (int*) malloc(lattice_size * sizeof(int));

    for (int j=0; j<lattice_size; j++){
        reticolo_aus[j] = false;
    }

    int nold = 0;
    int nnew = 1;
    int lc = 1;

    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());

    clusterx[0] = rx;
    clustery[0] = ry;
    reticolo_aus[pos(rx, ry, L)] = true;

    int nnx[q], nny[q];
    while (nold < nnew){
        for (int p = nold; p < nnew; p++){
            nearest(clusterx[p], clustery[p], nnx, nny, L);  //nn of cluster[p]
            for (int i = 0; i < q; i++){
                if (reticolo_aus[pos(nnx[i], nny[i], L)] == false && reticolo[pos(nnx[i], nny[i], L)] == reticolo[pos(clusterx[p], clustery[p], L)] && myrand() < prob){                    
                    clusterx[lc] = nnx[i];
                    clustery[lc] = nny[i];

                    reticolo_aus[pos(nnx[i], nny[i], L)] = true;
                    lc++;
                }
            }
        }

        nold = nnew;
        nnew = lc;
    }

    for (int i = 0; i < lc; i++){
        reticolo[pos(clusterx[i], clustery[i], L)] *= -1;
    }
    
    free(reticolo_aus);
    free(clusterx);
    free(clustery);

    return lc;
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

        fprintf(fp, "%f, %f\n", magn_ising(lattice, lattice_size), energy(lattice, lattice_size, L, beta));    //questo rallenta molto
        
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

    lattice_size = L*L;
    int *restrict lattice = (int *) malloc(lattice_size * sizeof(int));
    double prob = 1 -  exp(- 2 * beta);

    initialize_lattice_ising(lattice, lattice_size);
    
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
    int iter_bet_meas = 1;    //iterations between two measures
    int num_measures = 1e4;
    bool save_config = false;

    choose_geometry(modello, &nearest, &energy, &q);

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;
    myrand_init(seed1, seed2);
    
    
    int L_array[] = {40};
    double beta_array[] = {.67, 1};

    int num_L = sizeof(L_array) / sizeof(int);
    int num_beta = sizeof(beta_array) / sizeof(double);

    //omp_set_num_threads(2);
    //#pragma omp parallel for collapse(2)  //bisogna passare beta e L come input
    for (int i = 0; i < num_beta; i++){
        for (int j = 0; j < num_L; j++){
            beta = beta_array[i];
            L = L_array[j];

            montecarlo(iterations, iter_bet_meas, num_measures, save_config);
        }
    }

    return EXIT_SUCCESS;
}