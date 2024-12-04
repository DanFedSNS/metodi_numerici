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
/*void computejack(double *restrict datajack,
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
*/
void computejack(double *restrict datajack,
                 FILE *fp,
                 long int numberofbins,
                 int binsize,
                 double eta,
                 long int sampleeff,
                 int num_vars,
                 int skip_lines,
                 double g)
{
    long int i, r;
    int var;
    double *totals = malloc(num_vars * sizeof(double));   // To store totals for each variable
    double *subtotals = malloc(num_vars * sizeof(double)); // To store subtotals in the bin loop
    double *data_row = malloc(num_vars * sizeof(double)); // Temporary row to store a line of data

    // Initialize totals to zero
    for (var = 0; var < num_vars; var++) {
        totals[var] = 0.0;
    }

    char header[1000];
	rewind(fp); // Reset file pointer to re-read the relevant bin
    fgets(header, sizeof(header), fp); // Skip the header line
    // Skip the first 'skip_lines' lines
    for (i = 0; i < skip_lines; i++) {
        if (fgets(header, sizeof(header), fp) == NULL) {
            printf("\nIl valore di skip_lines = %d è maggiore della lunghezza del file\n", skip_lines);
            exit(1);
        }
    }
    // Calculate the total sums for all variables and track bin assignments
    for (i = 0; i < sampleeff; i++) {
        for (var = 0; var < num_vars; var++) {
            if (fscanf(fp, "%lf", &data_row[var]) != 1) {
                fprintf(stderr, "Error reading variable %d for data index %ld\n", var, i);
                exit(EXIT_FAILURE);
            }
            totals[var] += data_row[var];
        }        
    }

	rewind(fp); // Reset file pointer to re-read the relevant bin
    fgets(header, sizeof(header), fp); // Skip the header line
    // Skip the first 'skip_lines' lines
    for (i = 0; i < skip_lines; i++) {
        if (fgets(header, sizeof(header), fp) == NULL) {
            printf("\nIl valore di skip_lines = %d è maggiore della lunghezza del file\n", skip_lines);
            exit(1);
        }
    }
    // Loop over bins to compute jackknife subtotals
    for (i = 0; i < numberofbins; i++) {
        // Copy totals to subtotals
        for (var = 0; var < num_vars; var++) {
            subtotals[var] = totals[var];
        }

        // Subtract bin data from subtotals
        for (r = 0; r < binsize; r++) {
            for (var = 0; var < num_vars; var++) {
                if (fscanf(fp, "%lf", &data_row[var]) != 1) {
                    fprintf(stderr, "Error reading variable %d for bin index %ld\n", var, i);
                    exit(EXIT_FAILURE);
                }
                subtotals[var] -= data_row[var];
            }
        }

        // Normalize and assign values to `datajack`
        for (var = 0; var < num_vars; var++) {
            subtotals[var] /= (double)((numberofbins - 1) * binsize);
            datajack[(num_vars + 1) * i + var] = subtotals[var]; // Use an extended column width (n + 1)
        }

        // Compute the special value and place it in the last column
        if (num_vars >= 3) { // Ensure there are at least three variables
            datajack[(num_vars + 1) * i + num_vars] = 0.5 * subtotals[1] - subtotals[2] + 0.5 / eta + eta * g * subtotals[3];
        }
    }

    // Free allocated memory
    free(totals);
    free(subtotals);
    free(data_row);
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

void analysis(int Nt, double eta, double g){
	int skip_lines, binsize, j;
	long int sample, measevery, numberofbins, sampleeff, i;
	double *datajack;
	FILE *fp;
	char datafile[STRING_LENGTH];
	int num_vars;	//x, x^2, K, corrs
	skip_lines = 1e4;
    double simbeta = eta * (double)Nt;
    
	sprintf(datafile, "./misure/Nt%d_simbeta%.3f_g%.3f.dat", Nt, simbeta, g);
	num_vars = count_columns(datafile);
	
	fp = fopen(datafile, "r");
	char header[1000];
    fgets(header, sizeof(header), fp);  //ottiene la prima riga
	
    sscanf(header, "Nt = %*d, simbeta = %*f, g = %*f, sample = %ld, measevery = %ld, num_deltat = %*d", &sample, &measevery);	//attenzione, sample non è il numero di misure, bisogna dividere per measevery
    sample /= measevery;

	binsize = sample / 100;
	
	// definition of eta
	eta = simbeta / (double)Nt;
	
	// initialize numberofbins and sampleeff
	sample -= skip_lines;
	numberofbins = sample / binsize;
	sampleeff = numberofbins * binsize;
	
	// allocate jackknife samples
	datajack = (double *)malloc((unsigned long int)((num_vars + 1) * numberofbins) * sizeof(double)); // 4 because there are 4 observables
	if (datajack == NULL){
		fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
		exit(1);
	}
	
 
	// compute jackknife resamplings
    computejack(datajack, fp, numberofbins, binsize, eta, sampleeff, num_vars, skip_lines, g);
    fclose(fp);
    
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
	free(datajack);

	char datafile_o[50]; // file name
    sprintf(datafile_o, "./analysis/Nt%d_simbeta%.3f_g%.3f.dat", Nt, simbeta, g);
	
	fp = fopen(datafile_o, "w");
	fprintf(fp, "%s\n", header);

	for (j = 0; j < num_res; j++){
		fprintf(fp, "%.10f %.10f\n", ris[j], err[j]);
	}

	fclose(fp);
}

int main(){
    int Nt_start = 20;
    int Nt_stop = 100;
    int num_Nt = 1;
    int Nt_array[num_Nt];
    //arange_int(Nt_array, Nt_start, Nt_stop, num_Nt);
	
	
	double g_start = -2.0;
    double g_stop = 1.0;
    int num_g = 11;
    double g_array_aux[num_g];
    double g_array[num_g];
	
    linspace(g_array_aux, g_start, g_stop, num_g);
	for (int j = 0; j < num_g; j++){
        g_array[j] = pow(10, g_array_aux[j]);
    } 

	#pragma omp parallel for collapse(2) shared(Nt_array, g_array) schedule(dynamic, 1)  // collapse the loops and define private variables
	for (int i = 0; i < num_Nt; i++){
		for (int j = 0; j < num_g; j++){
			analysis(200, 1.0 / 28.0, g_array[j]);
		}
	}

	return EXIT_SUCCESS;
}