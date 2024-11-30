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

void update(int L, int *restrict reticolo, double *p, void (*nearest)(int, int, int *, int *, int), int q){
    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());

    int nnx[q], nny[q];
    nearest(rx, ry, nnx, nny, L);

    int Sr = 0;
    for (int i=0; i<q; i++){
        Sr += reticolo[pos(nnx[i], nny[i], L)];
    }

    int k = reticolo[pos(rx, ry, L)] * Sr;

    if (k <= 0){
        reticolo[pos(rx, ry, L)] *= -1;
    }

    else {
        double w = myrand();
        if (w <= p[k]){
            reticolo[pos(rx, ry, L)] *= -1;
        }
    }
}

void termalizzazione(int L, int *restrict lattice, int iterations, int lattice_size, double *p,
                     void (*nearest)(int, int, int *, int *, int), int q){
    //printf("\nAggiornamenti sulla termalizzazione:\n");
    for (int i=0; i < iterations * lattice_size; i++){  //termalizzazione
        update(L, lattice, p, nearest, q);

        /*if (i % (1000*lattice_size) == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)i / (double)lattice_size, (double)iterations);
        }*/
    }
    
}

void presa_misure(int L, double beta, int *restrict lattice, int lattice_size, int num_measures, 
                  int iter_bet_meas, FILE *fp, bool save_config, FILE *fp_config, double *p,
                  void (*nearest)(int, int, int *, int *, int), double (*energy)(int *restrict, int, int, double), int q){
    //printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas * lattice_size; i++){
            update(L, lattice, p, nearest, q);
        }

        fprintf(fp, "%.10f, %.10f\n", magn_ising(lattice, lattice_size), energy(lattice, lattice_size, L, beta));
        
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
    
    double p[7];     //pk = exp(-2beta k)
    for (int k=0; k <= q; k++){    //p[0] non verrà utilizzato mai, ma è comodo perché torna con k = sr * Sr
        p[k] = exp(- 2.0 * beta * ((double) k)); 
    }
    
    initialize_lattice_ising(lattice, lattice_size);
    
    FILE *fp, *fp_config; // pointer to file

    init_file(modello, L, beta, &fp, iterations, iter_bet_meas, num_measures, save_config, &fp_config);
    
    termalizzazione(L, lattice, iterations, lattice_size, p, nearest, q);
    
    presa_misure(L, beta, lattice, lattice_size, num_measures, iter_bet_meas, fp, save_config, fp_config, p, nearest, energy, q);
    
    free(lattice);

    close_file(&fp, save_config, &fp_config);

    double end = omp_get_wtime();
    double time_spent = (end - begin);
    //printf("\nHa impiegato %.3f secondi\n", time_spent);

    save_time_spent(beta, L, modello, time_spent, iterations, iter_bet_meas, num_measures);
}


int main(void){
    char *modello_values[] = {"ising_sq_metro"};//{"ising2d_tri_metro", "ising2d_sq_metro", "ising2d_hex_metro"};
    int num_modelli = sizeof(modello_values) / sizeof(modello_values[0]);

    int L_start = 10;
    int L_stop = 40;
    int num_L = 7;
    int L_array[num_L];
    arange_int(L_array, L_start, L_stop, num_L);

    int num_beta = 1;  
    
    for (int m = 0; m < num_modelli; m++) {
        #pragma omp parallel for collapse(2) shared(L_array, modello_values, num_beta) schedule(dynamic, 1) // collapse the loops and define private variables
        for (int i = 0; i < num_beta; i++) {
            for (int j = 0; j < num_L; j++) {
                int iterations = 0;
                int iter_bet_meas = 1;    //iterations between two measures
                int num_measures = 1e5;
                bool save_config = false;
                
                int q;  //nearest neighbours number
                void (*nearest)(int, int, int *, int *, int);
                double (*energy)(int *restrict, int, int, double);

                char *modello = modello_values[m];
                choose_geometry(modello, &nearest, &energy, &q);
                double beta = assign_beta(modello, i, num_beta);

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
