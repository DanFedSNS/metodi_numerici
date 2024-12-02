#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include"./include/get_array.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a two column file, each raw corresponding to measures taken on the same configuraton
// first column = x (averaged on Nt)
// second column = x^2 (averaged on Nt)
// thir columns = K (averaged on Nt)

// compute the jacknife samples of <x>, <x^2>, <Knaive>, <Hnaive>, <H>
void computejack(double *restrict datajack,
                 double const *const restrict data,
                 long int numberofbins,
                 int binsize,
                 double eta,
                 int n) // 'n' is the number of variables
{
    long int i, r;
    int j, var;
    const long int sampleeff = numberofbins * (long int)binsize;
    double *totals = malloc(n * sizeof(double)); // To store totals for each variable
    double *subtotals = malloc(n * sizeof(double)); // To store subtotals in the bin loop

    // Initialize totals to zero
    for (var = 0; var < n; var++) {
        totals[var] = 0.0;
    }

    // Calculate the total sums for all variables
    for (i = 0; i < sampleeff; i++) {
        for (var = 0; var < n; var++) {
            totals[var] += data[n * i + var];
        }
    }

    // Loop over bins
    for (i = 0; i < numberofbins; i++) {
        // Copy totals to subtotals
        for (var = 0; var < n; var++) {
            subtotals[var] = totals[var];
        }

        // Subtract bin data from subtotals
        for (j = 0; j < binsize; j++) {
            r = i * binsize + j;
            for (var = 0; var < n; var++) {
                subtotals[var] -= data[n * r + var];
            }
        }

        // Normalize and assign values to `datajack`
        for (var = 0; var < n; var++) {
            subtotals[var] /= (double)((numberofbins - 1) * binsize);
            datajack[(n + 1) * i + var] = subtotals[var]; // Use an extended column width (n + 1)
        }

        // Compute the special value and place it in the last column
        if (n >= 3) { // Ensure there are at least three variables
            datajack[(n + 1) * i + n] = 0.5 * subtotals[1] - subtotals[2] + 0.5 / eta;
        }
    }

    // Free allocated memory
    free(totals);
    free(subtotals);
}

int count_columns(char *datafile) {
	FILE *fp;
	fp = fopen(datafile, "r");
    char line[1024];
    int columns = 0;

	fgets(line, sizeof(line), fp);
	fgets(line, sizeof(line), fp);
    char *token;
	char *line_copy = strdup(line); // Create a copy of the line
	token = strtok(line_copy, " \t\n");
	while (token != NULL){
		columns++;
		token = strtok(NULL, " \t\n");
	}
	free(line_copy); // Free the duplicated line
	fclose(fp);
    return columns;
}

void analysis(int Nt, double simbeta){
	int skip_lines, binsize, j;
	long int sample, numberofbins, sampleeff, i;
	double eta;
	double *data, *datajack;
	FILE *fp;
	char datafile[STRING_LENGTH];
	int num_vars;	//x, x^2, K, corrs
	skip_lines = 1e4;

	sprintf(datafile, "./misure/Nt%d_simbeta%.1f.dat", Nt, simbeta);
	num_vars = count_columns(datafile);
	printf("num_vars = %d\n", num_vars);
	fp = fopen(datafile, "r");
	char header[1000];
    fgets(header, sizeof(header), fp);  //ottiene la prima riga
	
    sscanf(header, "Nt = %*d, simbeta = %*f, sample = %ld", &sample);	//attenzione, sample non è il numero di misure, bisogna dividere per measevery
 	
	binsize = sample / 100;
	
	// definition of eta
	eta = simbeta / (double)Nt;
	
	// initialize numberofbins and sampleeff
	sample -= skip_lines;
	numberofbins = sample / binsize;
	sampleeff = numberofbins * binsize;
	
	// allocate data arrays
	data = (double *)malloc((unsigned long int)(num_vars * sampleeff) * sizeof(double)); // 3 columns!
	if (data == NULL){
		fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
		exit(1);
	}

	// allocate jackknife samples
	datajack = (double *)malloc((unsigned long int)((num_vars + 1) * numberofbins) * sizeof(double)); // 4 because there are 4 observables
	if (datajack == NULL){
		fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
		exit(1);
	}
	
	// Skip the first 'skip_lines' lines
    for (i = 0; i < skip_lines; i++) {
        if (fgets(header, sizeof(header), fp) == NULL) {
            printf("\nIl valore di skip_lines = %d è maggiore della lunghezza del file\n", skip_lines);
            exit(1);
        }
    }

	for (i = 0; i < sampleeff; i++){
        for (int var = 0; var < num_vars; var++) {
			if (fscanf(fp, "%lf", &data[num_vars * i + var]) != 1){
				fprintf(stderr, "Error reading variable %d for data index %ld\n", var, i);
				exit(EXIT_FAILURE);
    		}
		}
	}


    fclose(fp);
     
	// compute jackknife resamplings
	computejack(datajack, data, numberofbins, binsize, eta, num_vars);

	int num_res = num_vars + 1;
	double ris[num_res], err[num_res];
	// compute average
	for (j = 0; j < num_res; j++){
		ris[j] = 0.0;
		for (i = 0; i < numberofbins; i++){
			ris[j] += datajack[num_res * i + j];
		}

		ris[j] /= (double)numberofbins;
	}

	// compute error
	for (j = 0; j < num_res; j++){
		err[j] = 0.0;
		for (i = 0; i < numberofbins; i++){
			err[j] += pow(ris[j] - datajack[num_res * i + j], 2.0);
		}
		// this corrects for a factor that is irrelevant but we leave it just for clarity
		err[j] *= (double)(numberofbins - 1);
		err[j] /= (double)numberofbins;
		err[j] = sqrt(err[j]);
	}

	// free data arrays
	free(data);
	free(datajack);

	char datafile_o[50]; // file name
    sprintf(datafile_o, "./analysis/Nt%d_simbeta%.1f.dat", Nt, simbeta);
	
	fp = fopen(datafile_o, "w");
	fprintf(fp, "Nt = %d, simbeta = %.10f\n", Nt, simbeta);

	for (j = 0; j < num_res; j++){
		fprintf(fp, "%.10f %.10f\n", ris[j], err[j]);
	}

	fclose(fp);
}

int main(){
	int L_start = 80;
    int L_stop = 190;
    int num_L = 12;
    int L_array[num_L];
    arange_int(L_array, L_start, L_stop, num_L);
	
	#pragma omp parallel for shared(L_array) schedule(dynamic, 1)  // collapse the loops and define private variables
	for (int i = 0; i < num_L; i++){
		analysis(L_array[i], 10.0);
	}

	

	return EXIT_SUCCESS;
}