#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "./include/random.h"
#include "./include/get_array.h"
#include <omp.h>

int L = 30;
const int D = 2;
double beta = .43;
int lattice_size;
double p[5];     //pk = exp(-2beta k)
int iterations = 2e5;
int iter_bet_meas = 1;    //iterations between two measures
int num_measures = 2e5;
bool save_config = false;

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

void update(int reticolo[L][L]){
    int rx = (int)((double)L * myrand());
    int ry = (int)((double)L * myrand());

    int nnx[2*D], nny[2*D];
    nearest(rx, ry, nnx, nny);

    int Sr = 0;
    for (int i=0; i<2*D; i++){
        Sr += reticolo[nnx[i]][nny[i]];
    }

    int k = reticolo[rx][ry] * Sr;

    if (k <= 0){
        reticolo[rx][ry] *= -1;
    }

    else {
        double w = myrand();
        if (w <= p[k]){
            reticolo[rx][ry] *= -1;
        }
    }

}

double magn(int reticolo[L][L]){
    int somma_spin = 0;
    for (int i=0; i < L;  i++){
        for (int j=0; j < L; j++){
            somma_spin += reticolo[i][j];
        }
    }
    
    return (double) somma_spin / (double) lattice_size;
}

double energy(int reticolo[L][L]){
  int tmp, sum;

  sum=0;
  for(int rx=0; rx<L; rx++){
     for(int ry=0; ry<L; ry++){
        tmp = (rx + 1) % L;
        sum += -reticolo[rx][ry] * reticolo[tmp][ry];
       
        tmp = (ry + 1) % L;
        sum += -reticolo[rx][ry] * reticolo[rx][tmp];
        }
     }

  return (double) sum / (double) (lattice_size);
}

void initialize_lattice(int lattice[L][L]){
    for (int i=0; i < L;  i++){
        for (int j=0; j < L; j++){
            lattice[i][j] = 1;
        } 
    }
}

void termalizzazione(int lattice[L][L], int iterations, int lattice_size){
    printf("\nAggiornamenti sulla termalizzazione:\n");
    for (int i=0; i < iterations * lattice_size; i++){  //termalizzazione
        update(lattice);

        if (i % (1000*lattice_size) == 0){
            printf("\33[2K\r");
            printf("Iterazione %.2e di %.2e", (double)i / (double)lattice_size, (double)iterations);
        }
    }
    
}

void presa_misure(int lattice[L][L], int lattice_size, int num_measures, int iter_bet_meas, FILE *fp, bool save_config, FILE *fp_config){
    printf("\nAggiornamenti sulle misure:\n");
    for (int j=0; j < num_measures; j++){    // presa misure
        for (int i=0; i < iter_bet_meas * lattice_size; i++){
            update(lattice);
        }

        fprintf(fp, "%f, %f\n", magn(lattice), energy(lattice));
        
        if (save_config == true){
            for (int a=0; a < L; a++){
                for (int b=0; b < L; b++){
                    fprintf(fp_config, "%d, ", lattice[a][b]);
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

void metropolis(){
    clock_t begin = clock();

    int lattice[L][L];
    lattice_size = L*L;

    for (int k=0; k <= 2*D; k++){    //p[0] non verrà utilizzato mai, ma è comodo perché torna con k = sr * Sr
        p[k] = exp(- 2.0 * beta * ((double) k)); 
    }
    
    initialize_lattice(lattice);
    
    char datafile[50], datafile_config[50]; // file name
    FILE *fp, *fp_config; // pointer to file

    sprintf(datafile, "./ising2d_metro/L%d_beta%.2f.dat", L, beta); // file name initialized with a string
    fp = fopen(datafile, "w");
    fprintf(fp, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);    
    
    if (save_config == true){
        sprintf(datafile_config, "./config/ising2d_metro_L%d_beta%.2f.dat", L, beta);
        fp_config = fopen(datafile_config, "w");
        fprintf(fp_config, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);
    }
    
    termalizzazione(lattice, iterations, lattice_size);
    
    presa_misure(lattice, lattice_size, num_measures, iter_bet_meas, fp, save_config, fp_config);
    
    fclose(fp);
    if (save_config == true){
        fclose(fp_config);    
    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nHa impiegato %.3f secondi\n", time_spent);


    char datafile_time[] = "./time/ising2d_metro";
    FILE *fp_time;

    fp_time = fopen(datafile_time, "a");
    fprintf(fp_time, "\ntime = %.2f min, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", time_spent/60, beta, L, iterations, iter_bet_meas, num_measures);    
    fclose(fp_time);
}


int main(void){
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;
    myrand_init(seed1, seed2);

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
    /*for (int i = 0; i < num_beta; i++){
        beta = beta_array[i];
        for (int j = 0; j < num_L; j++){
            L = L_array[j];

            metropolis();
        }
    }
    */
    
    return EXIT_SUCCESS;
}