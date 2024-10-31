#include <stdio.h>
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
