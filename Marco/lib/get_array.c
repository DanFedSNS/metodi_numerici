#include <stdio.h>
#include <stdbool.h>
#include "../include/get_array.h"

int count_lines(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error opening file");
        return -1; // Return -1 in case of an error
    }

    int line_count = 0; // Initialize line count
    char ch;

    // Read character by character until the end of the file
    while ((ch = fgetc(fp)) != EOF) {
        // Increment count for each newline character
        if (ch == '\n') {
            line_count++;
        }
    }

    fclose(fp); // Close the file
    return line_count; // Return the total line count
}

void get_array_from_txt_int(const char *datafile, int *array){
    FILE *fp;
    int num = count_lines(datafile);

    fp = fopen(datafile, "r");
    for (int i = 0; i < num; i++){
        fscanf(fp, "%d", &array[i]);
    }
    fclose(fp);
}

void get_array_from_txt_double(const char *datafile, double *array){
    FILE *fp;
    int num = count_lines(datafile);

    fp = fopen(datafile, "r");
    for (int i = 0; i < num; i++){
        fscanf(fp, "%lf", &array[i]);
    }
    fclose(fp);
}

void save_time_spent(double beta, int L, char *modello, double time_spent, int iterations, int iter_bet_meas, int num_measures){
    char datafile_time[50];
    sprintf(datafile_time, "./time/%s.dat", modello);
    FILE *fp_time;

    fp_time = fopen(datafile_time, "a");
    fprintf(fp_time, "\ntime = %.2f min, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", time_spent/60, beta, L, iterations, iter_bet_meas, num_measures);    
    fclose(fp_time);
}

void init_file(char *modello, int L, double beta, FILE **fp, int iterations, int iter_bet_meas, int num_measures, bool save_config, FILE **fp_config){
    char datafile[50]; // file name
    
    sprintf(datafile, "./%s/L%d_beta%.2f.dat", modello, L, beta); // file name initialized with a string
    *fp = fopen(datafile, "w");
    fprintf(*fp, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);    
    
    if (save_config == true){
        char datafile_config[50];
        sprintf(datafile_config, "./config/%s_L%d_beta%.2f.dat", modello, L, beta);
        *fp_config = fopen(datafile_config, "w");
        fprintf(*fp_config, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);
    }
} 

void close_file(FILE **fp, bool save_config, FILE **fp_config){
    fclose(*fp);
    if (save_config == true){
        fclose(*fp_config);    
    }
}

