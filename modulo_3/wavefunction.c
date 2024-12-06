#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include"./include/get_array.h"
#include "./include/boxmuller.h"
#include "./include/random.h"

#define STRING_LENGTH 100

// Metropolis update, return 1 if accepted
int metropolis(double *restrict lattice, long int r, double nnsum, double eta, const double delta, double g){
	double trial, Eold, Enew;

	Eold = lattice[r] * lattice[r] * (eta / 2.0 + 1. / eta) - lattice[r] * nnsum / eta + eta * g * pow(lattice[r], 4);
	trial = lattice[r] + delta * (1.0 - 2.0 * myrand());
	Enew = trial * trial * (eta / 2.0 + 1. / eta) - trial * nnsum / eta + eta * g * pow(trial, 4);

	if (Enew < Eold){
		lattice[r] = trial;
		return 1;
	}
	else if (myrand() < exp(-(Enew - Eold))){
		lattice[r] = trial;
		return 1;
	}

	return 0;
}

// overrelaxation
void overrelaxation(double *restrict lattice, long int r, double nnsum, double eta){
	const double average = nnsum / eta / (eta + 2.0 / eta);
	double ris = 2.0 * average - lattice[r];
	lattice[r] = ris;
}

int get_bin(double val, double xmin, double dx){
	return (int) ((val - xmin) / dx);
}

void save_time_spent(double time_spent, int Nt, double simbeta, double g, long int sample, int measevery, double acc_rate, int num_deltat){
	char datafile_time[50];
    sprintf(datafile_time, "./time/time.dat");
    FILE *fp_time;

    fp_time = fopen(datafile_time, "a");
    fprintf(fp_time, "\ntime = %.2f min, Nt = %d, simbeta = %.3f, g = %.3f, sample = %ld, measevery = %d, acc_rate = %f, num_deltat = %d\n",
						time_spent/60.0, Nt, simbeta, g, sample, measevery, acc_rate, num_deltat);    
    fclose(fp_time);	
}

int montecarlo(int Nt, double eta, long int sample, double g){
	double begin = omp_get_wtime();
	double *lattice;
	long int r, acc, bin_pos;
	long int *nnp, *nnm;
	double nnsum;
	double simbeta = eta * (double)Nt;
	double xmin = -3;
	double xmax = 3;
	int num_bins = 601;
	double dx = (xmax - xmin) / (double)(num_bins - 1); 
	long int hist[num_bins];

	for (r = 0; r < num_bins; r++){
		hist[r] = 0;
	}

	char datafile[STRING_LENGTH];
	sprintf(datafile, "./analysis/fig2/wf_Nt%d_simbeta%.3f_g%.3f.dat", Nt, simbeta, g); // file name initialized with a string

	FILE *fp;

	const int measevery = 100;
	//const int overrelaxsteps = 5;

	const double delta = 2.0 * sqrt(eta);

	// allocate the lattice and next neighbors
	lattice = (double *)malloc((unsigned long int)(Nt) * sizeof(double));
	if (lattice == NULL)
	{
		fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	nnp = (long int *)malloc((unsigned long int)(Nt) * sizeof(long int));
	if (nnp == NULL)
	{
		fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	nnm = (long int *)malloc((unsigned long int)(Nt) * sizeof(long int));
	if (nnm == NULL)
	{
		fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	// initialize nnp and nnm for periodic b.c.
	for (r = 0; r < Nt; r++)	{
		nnp[r] = r + 1;
		nnm[r] = r - 1;
	}
	nnp[Nt - 1] = 0;
	nnm[0] = Nt - 1;

	// initialize lattice
	for (r = 0; r < Nt; r++){
		lattice[r] = 0.0;
	}

	acc = 0;
	int out_bounds = 0;
	for (int iter = 0; iter < sample; iter++){
		for (int i = 0; i < Nt; i++){
			r = (int)((double)Nt * myrand());
			nnsum = lattice[nnp[r]] + lattice[nnm[r]];
			acc += metropolis(lattice, r, nnsum, eta, delta, g);
		}
	
		if (iter % measevery == measevery - 1){
			for (int i = 0; i < Nt; i++){
				bin_pos = get_bin(lattice[i], xmin, dx);
				if (bin_pos >= 0 && bin_pos < num_bins){
					hist[bin_pos] += 1;	
				} 
				else {
					out_bounds += 1;
				}
			}
		}
	}

	double acc_rate = (double)acc / (double)sample / (double)Nt;

	// open data file
	fp = fopen(datafile, "w");
	if (fp == NULL)	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	fprintf(fp, "Nt = %d, simbeta = %.10f, g = %.10f, sample = %ld, measevery = %d, out_bounds = %ld\n", Nt, simbeta, g, sample, measevery, out_bounds);

	for (r = 0; r < num_bins; r++){
		fprintf(fp, "%.10f %ld\n", xmin + r * dx, hist[r]);
	}
	
	fclose(fp);

	free(lattice);
	free(nnp);
	free(nnm);

	double end = omp_get_wtime();
    double time_spent = (end - begin);

	//save_time_spent(time_spent, Nt, simbeta, g, sample, measevery, acc_rate);

	return EXIT_SUCCESS;
}

int main(void){
	int Nt_start = 20;
    int Nt_stop = 100;
    int num_Nt = 1;
    int Nt_array[num_Nt];
    //arange_int(Nt_array, Nt_start, Nt_stop, num_Nt);
	
	
	double g_start = -2.0;
    double g_stop = 1.0;
    int num_g = 1;
    double g_array_aux[num_g];
    double g_array[num_g];
	
    linspace(g_array_aux, g_start, g_stop, num_g);
	for (int j = 0; j < num_g; j++){
        g_array[j] = pow(10, g_array_aux[j]);
    } 

	#pragma omp parallel for collapse(2) shared(Nt_array, g_array) schedule(dynamic, 1)  // collapse the loops and define private variables
	for (int i = 0; i < num_Nt; i++){
		for (int j = 0; j < num_g; j++){
			const unsigned long int seed1=(unsigned long int) time(NULL) + omp_get_thread_num();
        	const unsigned long int seed2=seed1+127;

        	// initialize random number generator
        	myrand_init(seed1, seed2);
			montecarlo(200, 1.0 / 28, 1e6, 0);
		}
	}

	return EXIT_SUCCESS;
}
