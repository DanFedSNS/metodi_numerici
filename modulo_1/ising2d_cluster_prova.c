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

const int D = 2;

void build_cluster_norec(int const * const restrict reticolo, 
                         long int r, 
                         int * restrict occup, 
                         long int * restrict pointtoocc, 
                         long int * restrict clustersize,
                         void (*nearest)(int, int, int *, int *, int),
                         long int volume,
                         double prob, int q, int L){
    long int r1, nnew;    
    long int nold=0; // starting value, with *clustersize=1

    // if first neighbors have the same orientation and are not occupied
    // they are added to the cluster with probability prob

    int nnx[q], nny[q];
    while(*clustersize>nold){ // this means that there are sites recently added, whose neighbors has not been checked yet, so we check them
        nnew=*clustersize;

        for (int p = nold; p < nnew; p++){
            r1=pointtoocc[p];
            nearest(r1/L, r1%L, nnx, nny, L);  //nn of cluster[p]
            
            for (int i = 0; i < q; i++){
                if(occup[pos(nnx[i], nny[i], L)]==0 && reticolo[r1]*reticolo[pos(nnx[i], nny[i], L)]==1){
                    if(myrand()<prob){
                        occup[pos(nnx[i], nny[i], L)]=1;
                        pointtoocc[*clustersize] = pos(nnx[i], nny[i], L);
                        (*clustersize)++;
                    }
                }
            }
        }

        nold=nnew;
    }
}

void update(int L, int *restrict reticolo, int lattice_size, double prob,
            void (*nearest)(int, int, int *, int *, int), int q){
    int *occup=(int *)malloc((unsigned long int)(lattice_size)*sizeof(int));
    long int *pointtoocc = (long int *)malloc((unsigned long int)(lattice_size)*sizeof(long int));
    
    for(int r=0; r<lattice_size; r++){
        occup[r]=0;
    }
    long int clustersize = 0;

    int r=(int)((double)lattice_size*myrand());
    occup[r]=1; // r is set as occupied
    pointtoocc[clustersize]=r; // a pointer to "r" is added in position "clustersize"
    clustersize++;

    build_cluster_norec(reticolo, r, occup, pointtoocc, &clustersize, nearest, lattice_size, prob, q, L);

    for(int r=0; r<clustersize; r++){
        reticolo[pointtoocc[r]] = -reticolo[pointtoocc[r]];
    }

    free(occup);
    free(pointtoocc);
}

void presa_misure(int L, double beta, int *restrict lattice, int lattice_size, int num_measures, int iter_bet_meas, 
                  FILE *restrict fp, FILE *restrict fp_config, bool save_config, double prob,
                  void (*nearest)(int, int, int *, int *, int), double (*energy)(int *restrict, int, int, double), int q){
    //printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas; i++){
            update(L, lattice, lattice_size, prob, nearest, q);
        }

        fprintf(fp, "%.10f, %.10f\n", magn_ising(lattice, lattice_size), energy(lattice, lattice_size, L, beta));    //questo rallenta molto
        
        if (save_config == true){
            for (int a=0; a < lattice_size; a++){
                fprintf(fp_config, "%d, ", lattice[a]);
            }       
            fseek(fp_config, -2, SEEK_CUR); //rimuove ", " finali
            fprintf(fp_config, "\n");
        }

        /*if (j % 1000 == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)j, (double)num_measures);
        }*/
    }
}

void montecarlo(int L, double beta, char *modello, int iterations, int iter_bet_meas, int num_measures, bool save_config,
                void (*nearest)(int, int, int *, int *, int), double (*energy)(int *restrict, int, int, double), int q){
    double begin = omp_get_wtime();

    int lattice_size = L*L;
    int *restrict lattice = (int *) malloc(lattice_size * sizeof(int));
    double prob = 1.0 - exp(- 2.0 * beta);

    initialize_lattice_ising(lattice, lattice_size);
    
    FILE *fp, *fp_config; //pointer to file

    init_file(modello, L, beta, &fp, iterations, iter_bet_meas, num_measures, save_config, &fp_config);

    presa_misure(L, beta, lattice, lattice_size, num_measures, iter_bet_meas, fp, fp_config, save_config, prob, nearest, energy, q);

    free(lattice);

    close_file(&fp, save_config, &fp_config);

    double end = omp_get_wtime();
    double time_spent = (end - begin);
    //printf("\nHa impiegato %.3f secondi\n", time_spent);

    save_time_spent(beta, L, modello, time_spent, iterations, iter_bet_meas, num_measures);
}


int main(void){
    char *modello_values[] = {"ising2d_sq_cluster"}; //{"ising2d_tri_cluster", "ising2d_sq_cluster", "ising2d_hex_cluster"};
    int num_modelli = sizeof(modello_values) / sizeof(modello_values[0]);
    
    int num_L = 1;
    int L_array[] = {90};

    int num_beta = 40;
    
    for (int m = 0; m < num_modelli; m++) {       
        #pragma omp parallel for collapse(2) shared(L_array, modello_values, num_beta) schedule(dynamic, 1)  // collapse the loops and define private variables
        for (int i = 0; i < num_beta; i++){
            for (int j = 0; j < num_L; j++){
                int iterations = 0;
                int iter_bet_meas = 1;    //iterations between two measures
                int num_measures = 1e6; 
                bool save_config = false; 

                int q;
                void (*nearest)(int, int, int *, int *, int);
                double (*energy)(int *restrict, int, int, double);

                char *modello = modello_values[m];
                choose_geometry(modello, &nearest, &energy, &q);
                double beta = assign_beta_close(modello, i, num_beta, L_array[j]);
                
                printf("\n%s, L = %d, beta = %.5f", modello, L_array[j], beta);
                // Initialize unique random seeds for each thread
                unsigned long int seed1 = (unsigned long int)time(NULL) + omp_get_thread_num();
                unsigned long int seed2 = seed1 + 127;
                myrand_init(seed1, seed2);

                montecarlo(L_array[j], beta, modello, iterations, iter_bet_meas, num_measures, save_config, nearest, energy, q);
            }
        }
    }

    return EXIT_SUCCESS;
}