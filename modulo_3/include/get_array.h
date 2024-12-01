#ifndef GET_ARRAY_H
#define GET_ARRAY_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// Function to count the number of lines in a file.
// Returns the number of lines or -1 if there's an error.
int count_lines(const char *filename);

// Function to read integers from a text file into an array.
// Assumes the array is large enough to hold all integers in the file.
void get_array_from_txt_int(const char *datafile, int *array);

// Function to read doubles from a text file into an array.
// Assumes the array is large enough to hold all doubles in the file.
void get_array_from_txt_double(const char *datafile, double *array);


void init_file(char *modello, int L, double beta, FILE **fp, int iterations, int iter_bet_meas, int num_measures, bool save_config, FILE **fp_config);

void close_file(FILE **fp, bool save_config, FILE **fp_config);


void nearest_sq(int rx, int ry, int *resx, int *resy, int L);  //nn reticolo quadrato

void nearest_tri(int rx, int ry, int *resx, int *resy, int L);  //nn reticolo triangolare

void nearest_hex(int rx, int ry, int *resx, int *resy, int L);  //nn reticolo esagonale (grafene)

void nearest_cu(int rx, int ry, int rz, int *resx, int *resy, int *resz, int L);  // nn reticolo cubico 

void initialize_lattice_ising(int *restrict lattice, int lattice_size);

double magn_ising(int *restrict reticolo, int lattice_size);

double energy_sq(int *restrict reticolo, int lattice_size, int L, double beta);

double energy_tri(int *restrict reticolo, int lattice_size, int L, double beta);

double energy_hex(int *restrict reticolo, int lattice_size, int L, double beta);

void choose_geometry(char *modello, void (**nearest)(int, int, int *, int *, int), double (**energy)(int *restrict, int, int, double), int *q);

void linspace(double *arr, double start, double stop, int num);

void arange_int(int *arr, int start, int stop, int num);

double assign_beta(char *modello, int i, int num_beta);

double assign_beta_close(char *modello, int i, int num_beta, int L);

#endif
