#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
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

void init_file(char *modello, int L, double beta, FILE **fp, int iterations, int iter_bet_meas, int num_measures, bool save_config, FILE **fp_config){
    char datafile[50]; // file name
    
    sprintf(datafile, "./%s/L%d_beta%.5f.dat", modello, L, beta); // file name initialized with a string
    *fp = fopen(datafile, "w");
    fprintf(*fp, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %d\n", beta, L, iterations, iter_bet_meas, num_measures);    
    
    if (save_config == true){
        char datafile_config[50];
        sprintf(datafile_config, "./config/%s_L%d_beta%.5f.dat", modello, L, beta);
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

void nearest_sq(int rx, int ry, int *resx, int *resy, int L){  //nn reticolo quadrato
    resx[0] = (rx - 1 + L) % L;     //+ L serve a evitare negativi
    resx[1] = (rx + 1) % L;
    resx[2] = rx;
    resx[3] = rx;

    resy[0] = ry;
    resy[1] = ry;
    resy[2] = (ry - 1 + L) % L;
    resy[3] = (ry + 1) % L;
}

void nearest_tri(int rx, int ry, int *resx, int *resy, int L){  //nn reticolo esagonale
    resx[0] = (rx - 1 + L) % L;     //+ L serve a evitare negativi
    resx[1] = (rx + 1) % L;
    resx[2] = rx;
    resx[3] = rx;
    resx[4] = (rx - 1 + L) % L;
    resx[5] = (rx + 1) % L;

    resy[0] = ry;
    resy[1] = ry;
    resy[2] = (ry - 1 + L) % L;
    resy[3] = (ry + 1) % L;
    resy[4] = (ry - 1 + L) % L;
    resy[5] = (ry + 1) % L;
}

void nearest_hex(int rx, int ry, int *resx, int *resy, int L){  //nn reticolo quadrato
    resx[0] = rx;     //+ L serve a evitare negativi
    resx[1] = rx;   //stessa cella
    if ((rx + ry) % 2 == 0){
        resx[2] = (rx + 1) % L;
    } else {
        resx[2] = (rx - 1 + L) % L;
    }

    resy[0] = (ry + 1) % L;
    resy[1] = (ry - 1 + L) % L;
    resy[2] = ry;
}

void nearest_cu(int rx, int ry, int rz, int *resx, int *resy, int *resz, int L){  // nn reticolo cubico 
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

void initialize_lattice_ising(int *restrict lattice, int lattice_size){
    for (int i=0; i < lattice_size;  i++){    
        lattice[i] = 1;
    }
}

double magn_ising(int *restrict reticolo, int lattice_size){
    int somma_spin = 0;
    for (int i=0; i < lattice_size;  i++){
            somma_spin += reticolo[i];
    }
    
    return (double) somma_spin / (double) lattice_size;
}

double energy_sq(int *restrict reticolo, int lattice_size, int L, double beta){
    #define pos_en(rx, ry, L) (rx * L + ry)     //ho cambiato nome per evitare problemi con pos
    double sum = 0.0;
    int rx_next, ry_next;
    
    for (int rx = 0; rx < L; rx++) {
        for (int ry = 0; ry < L; ry++) {
            // Interactions in the x-direction
            rx_next = (rx + 1) % L;
            sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx_next, ry, L)];

            // Interactions in the y-direction
            ry_next = (ry + 1) % L;
            sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx, ry_next, L)];
        }
    }

    return sum / (double) lattice_size;
} 

double energy_tri(int *restrict reticolo, int lattice_size, int L, double beta){
    #define pos_en(rx, ry, L) (rx * L + ry)     //ho cambiato nome per evitare problemi con pos
    double sum = 0.0;
    int rx_next, ry_next;
    
    for (int rx = 0; rx < L; rx++) {
        for (int ry = 0; ry < L; ry++) {
            // Interactions in the x-direction
            rx_next = (rx + 1) % L;
            sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx_next, ry, L)];

            // Interactions in the y-direction
            ry_next = (ry + 1) % L;
            sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx, ry_next, L)];

            sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx_next, ry_next, L)];
        }
    }

    return sum / (double) lattice_size;
} 

double energy_hex(int *restrict reticolo, int lattice_size, int L, double beta){
    #define pos_en(rx, ry, L) (rx * L + ry)     //ho cambiato nome per evitare problemi con pos
    double sum = 0.0;
    int rx_next, ry_next;
    
    for (int rx = 0; rx < L; rx++) {
        for (int ry = 0; ry < L; ry++) {        
            ry_next = (ry + 1) % L;
            sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx, ry_next, L)];

            if ((rx + ry) % 2 == 0){
                rx_next = (rx + 1) % L;
                sum += -reticolo[pos_en(rx, ry, L)] * reticolo[pos_en(rx_next, ry, L)];
            }
        }
    }

    return sum / (double) lattice_size;
}

void choose_geometry(char *modello, void (**nearest)(int, int, int *, int *, int), double (**energy)(int *restrict, int, int, double), int *q){
    if (strncmp(modello, "ising2d_sq", strlen("ising2d_sq")) == 0){
        *nearest = nearest_sq;
        *energy = energy_sq;
        *q = 4;
    }
    else if (strncmp(modello, "ising2d_tri", strlen("ising2d_tri")) == 0){
        *nearest = nearest_tri;
        *energy = energy_tri;
        *q = 6;
    }
    else if (strncmp(modello, "ising2d_hex", strlen("ising2d_hex")) == 0){
        *nearest = nearest_hex;
        *energy = energy_hex;
        *q = 3;
    }
    else {
        printf("\nIl nome del modello è sbagliato, il programma è stato interrotto.");
        exit(1);
    }
}

void linspace(double *arr, double start, double stop, int num){     //include start e stop, genera num elementi
    double delta = (stop - start) / (num - 1);
    for (int i = 0; i < num; i++){
        arr[i] = start + i * delta;
    }
}

void arange_int(int *arr, int start, int stop, int num){
    int delta = (stop - start) / (num - 1);
    for (int i = 0; i < num; i++){
        arr[i] = start + i * delta;
    }
}

double assign_beta(char *modello, int i, int num_beta){     //m indica il modello, i indica l'indice di beta_array
    double beta_start, beta_stop;

    if (strncmp(modello, "ising2d_tri", strlen("ising2d_tri")) == 0){
        beta_start = 0.27;
        beta_stop =  0.28;
    }
    else if (strncmp(modello, "ising2d_sq", strlen("ising2d_sq")) == 0){
        beta_start = 0.43;
        beta_stop = 0.45;        
    }
    else if (strncmp(modello, "ising2d_hex", strlen("ising2d_hex")) == 0){
        beta_start = 0.655;
        beta_stop = 0.685;
    }
    else{
        printf("\nErrore nella funzione assign_beta: modello = %s non valido\n", modello);
        exit(1);
    }

    double delta = (beta_stop - beta_start) / (num_beta - 1);
    return beta_start + i * delta;  //sarebbe l'i-esimo elemento di linspace(beta_start, beta_stop, num_beta)
}

double assign_beta_close(char *modello, int i, int num_beta, int L){   //close to beta_c
    double b_start, b_stop, beta_c;

    if (strncmp(modello, "ising2d_tri", strlen("ising2d_tri")) == 0){
        b_start = -0.5;
        b_stop =  0.125;
        beta_c = 0.2746531;
    }
    else if (strncmp(modello, "ising2d_sq", strlen("ising2d_sq")) == 0){
        b_start = -0.75;
        b_stop = 0.125;     
        beta_c = 0.4406868;   
    }
    else if (strncmp(modello, "ising2d_hex", strlen("ising2d_hex")) == 0){
        b_start = -1.0;
        b_stop = 0.5;
        beta_c = 0.6584789;
    }
    else{
        printf("\nErrore nella funzione assign_beta: modello = %s non valido\n", modello);
        exit(1);
    }

    double delta = (b_stop - b_start) / (num_beta - 1);
    double b_res = b_start + i * delta;  //sarebbe l'i-esimo elemento di linspace(beta_start, beta_stop, num_beta)

    return b_res / L + beta_c;
}

