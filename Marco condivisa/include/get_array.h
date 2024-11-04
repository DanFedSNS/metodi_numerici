#ifndef GET_ARRAY_H
#define GET_ARRAY_H

#include <stdio.h>
#include <stdbool.h>

// Function to count the number of lines in a file.
// Returns the number of lines or -1 if there's an error.
int count_lines(const char *filename);

// Function to read integers from a text file into an array.
// Assumes the array is large enough to hold all integers in the file.
void get_array_from_txt_int(const char *datafile, int *array);

// Function to read doubles from a text file into an array.
// Assumes the array is large enough to hold all doubles in the file.
void get_array_from_txt_double(const char *datafile, double *array);

void save_time_spent(double beta, int L, char *modello, double time_spent, int iterations, int iter_bet_meas, int num_measures);

void init_file(char *modello, int L, double beta, FILE **fp, int iterations, int iter_bet_meas, int num_measures, bool save_config, FILE **fp_config);

void close_file(FILE **fp, bool save_config, FILE **fp_config);


#endif
