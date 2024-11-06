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
int q;  //nearest neighbours number
double beta;
int lattice_size;
double p[7];     //pk = exp(-2beta k)
int iterations = 1e5;
int iter_bet_meas = 1;    //iterations between two measures
int num_measures = 1e4;
bool save_config = false;
void (*nearest)(int, int, int *, int *, int);
double (*energy)(int *restrict, int, int, double);
char modello[] = "ising2d_tri_metro";   //ising2d_tri_metro or ising2d_metro

void update(int *restrict reticolo){
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

void termalizzazione(int *restrict lattice, int iterations, int lattice_size){
    printf("\nAggiornamenti sulla termalizzazione:\n");
    for (int i=0; i < iterations * lattice_size; i++){  //termalizzazione
        update(lattice);

        if (i % (1000*lattice_size) == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)i / (double)lattice_size, (double)iterations);
        }
    }
    
}

void presa_misure(int *restrict lattice, int lattice_size, int num_measures, int iter_bet_meas, FILE *fp, bool save_config, FILE *fp_config){
    printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas * lattice_size; i++){
            update(lattice);
        }

        fprintf(fp, "%f, %f\n", magn_ising(lattice, lattice_size), energy(lattice, lattice_size, L, beta));
        
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

void montecarlo(){
    clock_t begin = clock();

    lattice_size = L*L;
    int *restrict lattice = (int *) malloc(lattice_size * sizeof(int));
    
    for (int k=0; k <= q; k++){    //p[0] non verrà utilizzato mai, ma è comodo perché torna con k = sr * Sr
        p[k] = exp(- 2.0 * beta * ((double) k)); 
    }
    
    initialize_lattice_ising(lattice, lattice_size);
    
    FILE *fp, *fp_config; // pointer to file

    init_file(modello, L, beta, &fp, iterations, iter_bet_meas, num_measures, save_config, &fp_config);
    
    termalizzazione(lattice, iterations, lattice_size);
    
    presa_misure(lattice, lattice_size, num_measures, iter_bet_meas, fp, save_config, fp_config);
    
    free(lattice);

    close_file(&fp, save_config, &fp_config);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nHa impiegato %.3f secondi\n", time_spent);

    save_time_spent(beta, L, modello, time_spent, iterations, iter_bet_meas, num_measures);
}


int main(void){
    if (strcmp(modello, "ising2d_metro") == 0){
        nearest = nearest_sq;
        energy = energy_sq;
        q = 4;
    }
    else if (strcmp(modello, "ising2d_tri_metro") == 0){
        nearest = nearest_tri;
        energy = energy_tri;
        q = 6;
    } else {
        printf("\nIl nome del modello è sbagliato, il programma è stato interrotto.");
        exit(1);
    }

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
    */

    int L_array[] = {40};
    double beta_array[] = {.4};

    int num_L = sizeof(L_array) / sizeof(int);
    int num_beta = sizeof(beta_array) / sizeof(double);

    //omp_set_num_threads(2);
    //#pragma omp parallel for collapse(2)  bisogna passare  beta e L come input
    for (int i = 0; i < num_beta; i++){
        for (int j = 0; j < num_L; j++){
            beta = beta_array[i];
            L = L_array[j];

            montecarlo();
        }
    }
    
    return EXIT_SUCCESS;
}