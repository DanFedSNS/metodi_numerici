#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<time.h>
#include "./include/get_array.h"
//#include<omp.h>
#define scarto(x, y) ((y - x) / x)

const int D = 2;

double autocorr(double *x, double x_avg, double x_std, int n, int len_x){   //funzione di autocorrelazione
    double res = 0;
    for (int i=0; i < len_x - n; i++){
        res += (x[i] - x_avg) * (x[i+n] - x_avg);
    }
    res *= 1.0 / pow(x_std, 2) / (len_x);
    return res;
}

double average(double *x, int len_x){   // <x>  -->  calcola la media di x
    double res = 0;
    for (int i = 0; i < len_x; i++){
        res += x[i];
    }

    return res / len_x;
}

double var(double *x, double x_avg, int len_x){  // <x^2> - <x>^2   --> varianza scorrelata
    double res = 0;
    for (int i = 0; i < len_x; i++){
        res += pow(x[i] - x_avg, 2);
    }

    return res / (len_x - 1);
}

double sd(double *x, double x_avg, int len_x){  // sqrt(<x^2> - <x>^2)   --> deviazione standard scorrelata
    double res = 0;
    for (int i = 0; i < len_x; i++){
        res += pow(x[i] - x_avg, 2);
    }

    return sqrt(res / (len_x - 1));
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

double binning(double *x, int len_x, int bin_min, int bin_max){   
    double x_avg = average(x, len_x);    // la calcolo una sola volta
    double *sigma_vals = (double *) malloc((bin_max + 1) * sizeof(double));
    double tolleranza_binning = .1;

    for (int k = bin_min; k <= bin_max; k++){
        int len_val_bin = len_x / k;    //quanti bin lunghi k ci sono
        double * restrict valori_nei_bin = malloc(len_val_bin * sizeof(double));
        raggruppa(x, valori_nei_bin, len_x, k);
        
        sigma_vals[k] = sd(valori_nei_bin, x_avg, len_val_bin) / sqrt(len_val_bin);
        
        free(valori_nei_bin);
    }

    double sigma_fin = 0.0; //sigma extracted from binning
    int values_considered;
    for (int k = bin_max; k > bin_min; k--){
        if (scarto(sigma_vals[k], sigma_vals[bin_max]) > tolleranza_binning){
            values_considered = bin_max - k;
            //printf("\nstopped at k = %d, bin_max = %d\n", k, bin_max);
            break;
        }
        
        sigma_fin += sigma_vals[k];
    }

    free(sigma_vals);

    return sigma_fin / values_considered;
}

// compute the jacknife samples of <E>, <E^2>-<E>^2, <M>, <|M|>, <M^2>-<|M|^2>, <M2>, <M^4>/<M^2>^2
void jacknife(double * restrict datajack, 
                 double const * const restrict magn, 
                 double const * const restrict energy, 
                 long int num_bins, 
                 int binsize){
    long int i, r;
    int j;
    const long int sampleeff=num_bins*(long int) binsize;
    double Etot, E2tot, Mabstot, M2tot;
    double E, E2, Mabs, M2;

    Etot=0.0;
    E2tot=0.0;
    Mabstot=0.0;
    M2tot=0.0;

    for(i=0; i<sampleeff; i++){
        Etot+=energy[i];
        E2tot+=pow(energy[i], 2.0);
        Mabstot+=magn[i];
        M2tot+=pow(magn[i],2.0);
    }

    for(i=0; i<num_bins; i++){
        E=Etot;
        E2=E2tot;
        Mabs=Mabstot;
        M2=M2tot;

        for(j=0; j<binsize; j++){
            r=i*binsize+j;

            E-=energy[r];
            E2-=pow(energy[r], 2.0);
            Mabs-=magn[r];
            M2-=pow(magn[r],2.0);
        }

        E/=(double)((num_bins-1)*binsize); 
        E2/=(double)((num_bins-1)*binsize);
        Mabs/=(double)((num_bins-1)*binsize);
        M2/=(double)((num_bins-1)*binsize);

        datajack[4*i+0]=E;
        datajack[4*i+1]=E2-E*E;
        datajack[4*i+3]=Mabs;
        datajack[4*i+4]=M2-Mabs*Mabs;
    }
}

void clear_initial_data(const char* datafile_o) {
    
    FILE *fp = fopen(datafile_o, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile_o, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    fclose(fp);
}

void analysis(int L, double beta, int skip_lines, char *modello){
    clock_t begin = clock();

    int lattice_size = (int) pow(L, D);
    char datafile[50]; 
    FILE *fp; // pointer to file
    sprintf(datafile, "./%s/L%d_beta%.5f.dat", modello, L, beta);
    fp = fopen(datafile, "r");

    if(fp==NULL){
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__); 
        exit(EXIT_FAILURE);
    }

    char header[1000];
    fgets(header, sizeof(header), fp);  //ottiene la prima riga
    
    int iterations, iter_bet_meas, num_measures;
    sscanf(header, "m, E, beta = %lf, L = %*d, iterations = %d, iter_bet_meas = %d, num_measures = %d", &beta, &iterations, &iter_bet_meas, &num_measures);
    
    // Skip the first 'skip_lines' lines
    for (int i = 0; i < skip_lines; i++) {
        if (fgets(header, sizeof(header), fp) == NULL) {
            printf("\nIl valore di skip_lines = %d è maggiore della lunghezza del file\n", skip_lines);
            exit(1);
        }
    }

    num_measures -= skip_lines;

    double * restrict magn = malloc(num_measures * sizeof(double));
    double * restrict energy = malloc(num_measures * sizeof(double));
    if (magn == NULL || energy == NULL) {   //check if malloc worked
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE); 
    }

    for (int i = 0; i < num_measures; i++){
        fscanf(fp, "%lf, %lf", &magn[i], &energy[i]);
        magn[i] = fabs(magn[i]);
    }

    fclose(fp);
    
    double energy_avg = average(energy, num_measures);
    double magn_abs_avg = average(magn, num_measures);
    double binder_cum = U(magn, num_measures);

    double specific_heat = pow(beta, 2) * lattice_size * var(energy, energy_avg, num_measures);
    double susceptibility = beta * lattice_size * var(magn, magn_abs_avg, num_measures);
    
    double * restrict chi_arr = malloc(num_measures * sizeof(double));
    double * restrict sp_heat_arr = malloc(num_measures * sizeof(double));

    for (int i = 0; i < num_measures; i++){
        chi_arr[i] = beta * lattice_size * pow(magn[i] - magn_abs_avg, 2);
        sp_heat_arr[i] = pow(beta, 2) * lattice_size * pow(energy[i] - energy_avg, 2);
    }
    
    int bin_max = num_measures / 50;
    int bin_min = 3;
    double sigma_magn = binning(magn, num_measures, bin_min, bin_max);
    double sigma_energy = binning(energy, num_measures, bin_min, bin_max);
    double sigma_susceptibility = binning(chi_arr, num_measures, bin_min, bin_max);
    double sigma_sp_heat = binning(sp_heat_arr, num_measures, bin_min, bin_max);
    
    free(chi_arr);
    free(sp_heat_arr);

    char datafile_o[50]; // file name
    sprintf(datafile_o, "./analysis_%s/L%d.dat", modello, L);

    fp = fopen(datafile_o, "a");

    fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
        beta, specific_heat, susceptibility, magn_abs_avg, energy_avg,
        binder_cum, sigma_magn, sigma_energy, sigma_susceptibility, sigma_sp_heat);

    
    fclose(fp);

    /*char datafile_ac[50]; // file name
    sprintf(datafile_ac, "./analysis_%s/L%d_autocorr.dat", modello, L);
    fp = fopen(datafile_ac, "a");

    //Autocorrelation function
    double energy_sd = sd(energy, energy_avg, num_measures);
    double magn_sd = sd(magn, magn_abs_avg, num_measures);
    int max_autocorr = 1e4;
    fprintf(fp, "autocorrelation_energy = ");
    for (int n = 0; n < max_autocorr; n++){
        fprintf(fp, "%lf, ", autocorr(energy, energy_avg, energy_sd, n, num_measures));
    }
    fseek(fp, -2, SEEK_CUR); //rimuove ", " finali
    
    for (int n = 0; n < max_autocorr; n++) {
        if (n == max_autocorr - 1) {
            fprintf(fp, "%lf", autocorr(magn, magn_abs_avg, magn_sd, n, num_measures));
        } else {
            fprintf(fp, "%lf, ", autocorr(magn, magn_abs_avg, magn_sd, n, num_measures));
        }
    }
    fprintf(fp, "\n");

    
    fclose(fp);*/
    free(magn);
    free(energy);
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("\nL = %d, beta = %f ha impiegato %.3f secondi\n", L, beta, time_spent);
}

int main(void) {
    char *modello_values[] = {"ising2d_tri_cluster"}; //, "ising2d_sq_cluster", "ising2d_hex_cluster"};
    int num_modelli = sizeof(modello_values) / sizeof(modello_values[0]);
    
    int L_start = 30;
    int L_stop = 40;
    int num_L = 2;
    int L_array[num_L];
    arange_int(L_array, L_start, L_stop, num_L);

    int num_beta = 2;
    int skip_lines = 0;
   
    #pragma omp parallel for collapse(2) shared(L_array, modello_values, num_beta, skip_lines)  // collapse the loops and define private variables
    for (int m = 0; m < num_modelli; m++){
        for (int j = 0; j < num_L; j++) {
            char *modello = modello_values[m];
            char datafile_o[50];
            char datafile_ac[50];
            
            sprintf(datafile_o, "./analysis_%s/L%d.dat", modello, L_array[j]); 
            clear_initial_data(datafile_o);

            sprintf(datafile_ac, "./analysis_%s/L%d_autocorr.dat", modello, L_array[j]);
            clear_initial_data(datafile_ac);

            for (int i = 0; i < num_beta; i++) {
                double beta = assign_beta_close(modello, i, num_beta, L_array[j]);
                analysis(L_array[j], beta, skip_lines, modello);
            }
        }
    }

    return EXIT_SUCCESS;
}