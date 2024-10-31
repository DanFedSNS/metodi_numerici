#ifndef GET_ARRAY_H
#define GET_ARRAY_H

#include <stdio.h>

// Function to count the number of lines in a file.
// Returns the number of lines or -1 if there's an error.
int count_lines(const char *filename);

// Function to read integers from a text file into an array.
// Assumes the array is large enough to hold all integers in the file.
void get_array_from_txt_int(const char *datafile, int *array);

// Function to read doubles from a text file into an array.
// Assumes the array is large enough to hold all doubles in the file.
void get_array_from_txt_double(const char *datafile, double *array);

#endif
