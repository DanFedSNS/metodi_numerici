#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<time.h>
#include "./include/get_array.h"
#include<omp.h>

int D = 3;
char modello[] = "ising2d_metro";

double autocorr(double *x, double x_avg, double x_std, int n, int len_x){   //funzione di autocorrelazione
    double res = 0;
    for (int i=0; i < len_x - n; i++){
        res += (x[i] - x_avg) * (x[i+n] - x_avg);
    }
    res *= 1.0 / pow(x_std, 2) / (len_x);
    return res;
}

double average(double *x, int len_x){   //calcola la media di x
    double res = 0;
    for (int i = 0; i < len_x; i++){
        res += x[i];
    }

    return res / len_x;
}

double abs_average(double *x, int len_x){   //<|x|>
    double res = 0;
    for (int i = 0; i < len_x; i++){
        res += fabs(x[i]);
    }

    return res / len_x;
}

double sd(double *x, double x_avg, int len_x){  //deviazione standard di x scorrelata
    double res = 0;
    for (int i = 0; i < len_x; i++){
        res += pow(x[i] - x_avg, 2);
    }

    return sqrt(res / (len_x - 1));
}

double chi(double *x, int len_x){  //<x^2> - <|x|>^2
    double res = 0;
    double abs_avg = 0;

    for (int i = 0; i < len_x; i++){
        res += pow(x[i], 2);
        abs_avg += fabs(x[i]);
    }

    res /= len_x;
    abs_avg /= len_x;
    
    return res - pow(abs_avg, 2);
}

double U(double *x, int len_x){
    double numeratore = 0;
    double denominatore = 0;

    for (int i = 0; i < len_x; i++){
        numeratore += pow(x[i], 4);
        denominatore += pow(x[i], 2);
    }

    return numeratore / pow(denominatore, 2) * len_x;
}

void raggruppa(double *x, double *binned_x, int len_x, int k_bin){  //prende x e lo divide in bin lunghi k
    for (int i = 0; i < len_x / k_bin; i++){
        binned_x[i] = 0;
        for (int j = 0; j < k_bin; j++){
            binned_x[i] += x[k_bin * i + j];
        }
        binned_x[i] /= k_bin;
    }
}

void binning(double *x, double *sigma_vals, int len_x, int bin_min, int bin_max){     //sigma_vals conterrà i valori di std(F(k))
    double x_avg = average(x, len_x);    // la calcolo una sola volta

    for (int k = bin_min; k <= bin_max; k++){
        int len_val_bin = len_x / k;    //quanti bin lunghi k ci sono
        double * restrict valori_nei_bin = malloc(len_val_bin * sizeof(double));
        raggruppa(x, valori_nei_bin, len_x, k);
        
        sigma_vals[k] = sd(valori_nei_bin, x_avg, len_val_bin) / sqrt(len_val_bin);
        
        free(valori_nei_bin);
    }
}

void analysis(int L, double beta){
    clock_t begin = clock();
    
    char datafile[50]; // file name
    FILE *fp; // pointer to file
    sprintf(datafile, "./%s/L%d_beta%.2f.dat", modello, L, beta);
    fp = fopen(datafile, "r");

    if(fp==NULL){
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__); 
        exit(EXIT_FAILURE);
    }

    char header[1000];
    fgets(header, sizeof(header), fp);  //ottiene la prima riga
    
    int iterations, iter_bet_meas, num_measures;
    sscanf(header, "m, E, beta = %lf, L = %*d, iterations = %d, iter_bet_meas = %d, num_measures = %d", &beta, &iterations, &iter_bet_meas, &num_measures);
    
    double * restrict magn = malloc(num_measures * sizeof(double));
    double * restrict energy = malloc(num_measures * sizeof(double));
    if (magn == NULL || energy == NULL) {   //check if malloc worked
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE); 
    }

    
    for (int i = 0; i < num_measures; i++){
        fscanf(fp, "%lf, %lf", &magn[i], &energy[i]);
    }

    fclose(fp);
    
    double energy_avg = average(energy, num_measures);
    double energy_sd = sd(energy, energy_avg, num_measures);
    double magn_abs_avg = abs_average(magn, num_measures);
    double magn_sd = sd(magn, magn_abs_avg, num_measures);
    double binder_cum = U(magn, num_measures);

    double specific_heat = pow(beta, 2) * pow(L,D) * pow(energy_sd, 2);
    double susceptibility = beta * pow(L,D) * chi(magn, num_measures);
    /*
    //binning-blocking
    int bin_max = num_measures / 50;
    int bin_min = 3;
    int bin_num = bin_max - bin_min + 1;
    double * restrict sigma_magn_binned = malloc(bin_num * sizeof(double));
    double * restrict sigma_energy_binned = malloc(bin_num * sizeof(double));
    
    binning(magn, sigma_magn_binned, num_measures, bin_min, bin_max);
    binning(energy, sigma_energy_binned, num_measures, bin_min, bin_max);
    */

    //output data to file
    char datafile_o[50]; // file name
    sprintf(datafile_o, "./analysis_%s/L%d_beta%.2f.dat", modello, L, beta); // file name initialized with a string

    fp = fopen(datafile_o, "w");
    
    fprintf(fp, "beta = %lf, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);
    fprintf(fp, "specific_heat = %lf\n", specific_heat);
    fprintf(fp, "susceptibility = %lf\n", susceptibility);
    fprintf(fp, "magn_abs_avg = %lf\n", magn_abs_avg);
    fprintf(fp, "energy_avg = %lf\n", energy_avg);
    fprintf(fp, "binder_cum = %lf\n", binder_cum);
    
    /*
    fprintf(fp, "sigma_energy_bin = ");
    for (int i = 0; i < bin_num; i++){
        fprintf(fp, "%lf, ", sigma_energy_binned[i]);
    }
    fseek(fp, -2, SEEK_CUR); //rimuove ", " finali
    
    
    fprintf(fp, "\nsigma_magn_bin = ");
    for (int i = 0; i < bin_num; i++){
        fprintf(fp, "%lf, ", sigma_magn_binned[i]);
    }
    fseek(fp, -2, SEEK_CUR); //rimuove ", " finali
    
    free(sigma_energy_binned);
    free(sigma_magn_binned);
    
    //Autocorrelation function  INUTILE SCRIVERE TUTTI QUESTI DATI, NE BASTA 1 OGNI e.g. 10
    int max_autocorr = num_measures / 10;
    fprintf(fp, "\nautocorrelation_energy = ");
    for (int n = 0; n < max_autocorr; n++){
        fprintf(fp, "%lf, ", autocorr(energy, energy_avg, energy_sd, n, num_measures));
    }
    fseek(fp, -2, SEEK_CUR); //rimuove ", " finali
    
    fprintf(fp, "\nautocorrelation_magn = ");
    for (int n = 0; n < max_autocorr; n++){
        fprintf(fp, "%lf, ", autocorr(magn, magn_abs_avg, magn_sd, n, num_measures));
    }
    fseek(fp, -2, SEEK_CUR); //rimuove ", " finali
    */

    fprintf(fp, "  ");  //la riga sopra sposta solo il cursore, se non sovrascrivi ", " rimane
    
    fclose(fp);
    free(magn);
    free(energy);
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nL = %d, beta = %f ha impiegato %.3f secondi\n", L, beta, time_spent);
}


int main(void){
    /*
    const char *datafile_Larray, *datafile_betaarray;
    datafile_Larray = "L_array.txt";
    datafile_betaarray = "beta_array.txt";

    int num_L = count_lines(datafile_Larray);
    int num_beta = count_lines(datafile_betaarray);

    int L_array[num_L];
    double beta_array[num_beta];

    get_array_from_txt_int(datafile_Larray, L_array);
    get_array_from_txt_double(datafile_betaarray, beta_array);
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < num_beta; i++){
        for (int j = 0; j < num_L; j++){
            printf("i, j = %d, %d\n", i, j);
            analysis(L_array[j], beta_array[i]);
        }
    }
    */
    
    int L = 30;
    double beta = 0.23;
    analysis(L, beta);
    
    return EXIT_SUCCESS;

}