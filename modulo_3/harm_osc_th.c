#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "./include/boxmuller.h"
#include "./include/random.h"

#define STRING_LENGTH 50

#define METRO // if this macro is defined Metropolis is uded, otherwise heatbath

// average position
double calc_x(double const *const restrict lattice, long int Nt)
{
	long int r;
	double ris;

	ris = 0.0;
	for (r = 0; r < Nt; r++)
	{
		ris += lattice[r];
	}

	return ris / (double)Nt;
}

// average position square
double calc_x2(double const *const restrict lattice, long int Nt)
{
	long int r;
	double ris;

	ris = 0.0;
	for (r = 0; r < Nt; r++)
	{
		ris += lattice[r] * lattice[r];
	}

	return ris / (double)Nt;
}

// <x^n>
double calc_xn(double const *const restrict lattice, long int Nt, int n)
{
	long int r;
	double res = 0.0;

	for (r = 0; r < Nt; r++)
	{
		res += pow(lattice[r], n);
	}

	return res / (double)Nt;
}

// average kinetic energy (naive!)
double calc_Knaive(double const *const restrict lattice,
				   long int const *const restrict nnp,
				   long int Nt,
				   double eta)
{
	long int r;
	double ris, aux;

	ris = 0.0;
	for (r = 0; r < Nt; r++)
	{
		aux = lattice[nnp[r]] - lattice[r];
		ris += aux * aux / eta / eta;
	}

	return ris / (double)Nt / 2.0;
}

// <x^n(deltat) x^n(0)> correlator
double correlator(double const *const restrict lattice, long int Nt, long int deltat, int n, int m)
{
	long int r, raux;
	int i;
	double res = 0;

	for (r = 0; r < Nt; r++){
		raux = (r + deltat) % Nt;

		res += pow(lattice[r], n) * pow(lattice[raux], m);
	}

	return res / (double)Nt;
}

// Metropolis update, return 1 if accepted
int metropolis(double *restrict lattice, long int r, double nnsum, double eta, const double delta)
{
	double trial, Eold, Enew;

	Eold = lattice[r] * lattice[r] * (eta / 2.0 + 1. / eta) - lattice[r] * nnsum / eta;
	trial = lattice[r] + delta * (1.0 - 2.0 * myrand());
	Enew = trial * trial * (eta / 2.0 + 1. / eta) - trial * nnsum / eta;

	if (Enew < Eold)
	{
		lattice[r] = trial;
		return 1;
	}
	else if (myrand() < exp(-(Enew - Eold)))
	{
		lattice[r] = trial;
		return 1;
	}

	return 0;
}

// Heatbath update, return 1
int heatbath(double *restrict lattice, long int r, double nnsum, double eta)
{
	const double std = 1.0 / sqrt(eta + 2.0 / eta);
	const double average = nnsum / eta / (eta + 2.0 / eta);

	lattice[r] = average + std * gauss1();

	return 1;
}

// overrelaxation
void overrelaxation(double *restrict lattice, long int r, double nnsum, double eta)
{
	const double average = nnsum / eta / (eta + 2.0 / eta);
	double ris = 2.0 * average - lattice[r];
	lattice[r] = ris;
}

int montecarlo(int Nt, double simbeta, long int sample)
{
	double begin = omp_get_wtime();
	double *lattice;
	long int r, acc;
	long int *nnp, *nnm;
	double nnsum;
	double x, x2, Knaive;

	char datafile[] = "prova.dat";
	FILE *fp;

	const int measevery = 10;
	const int overrelaxsteps = 1;

	double eta = simbeta / (double)Nt;
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
	for (r = 0; r < Nt; r++)
	{
		nnp[r] = r + 1;
		nnm[r] = r - 1;
	}
	nnp[Nt - 1] = 0;
	nnm[0] = Nt - 1;

	// initialize lattice
	for (r = 0; r < Nt; r++)
	{
		lattice[r] = 0.0;
	}

	// open data file
	fp = fopen(datafile, "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	acc = 0;
	for (int iter = 0; iter < sample; iter++)
	{
		if (myrand() < 0.5)
		{
			for (r = 0; r < Nt; r++)
			{
				nnsum = lattice[nnp[r]] + lattice[nnm[r]];

				#ifdef METRO
					acc += metropolis(lattice, r, nnsum, eta, delta);
				#else
				acc += heatbath(lattice, r, nnsum, eta);
				#endif
			}
		}
		else
		{
			// overrelaxation
			for (int j = 0; j < overrelaxsteps; j++){
				for (r = 0; r < Nt; r++){
					nnsum = lattice[nnp[r]] + lattice[nnm[r]];

					overrelaxation(lattice, r, nnsum, eta);
				}
			}
		}

		if (iter % measevery == 0){
			x = calc_x(lattice, Nt);
			x2 = calc_x2(lattice, Nt);
			Knaive = calc_Knaive(lattice, nnp, Nt, eta);

			fprintf(fp, "%.10f %.10f %.10f ", x, x2, Knaive);
			int n_corr_max = 2;
			for (int n = 1; n <= n_corr_max; n++){
				for (int m = 1; m <= n_corr_max; m++){
					for (int deltat = 0; deltat <= Nt / 4; deltat+= Nt/4){
						fprintf(fp, "%.10f ", correlator(lattice, Nt, deltat, n, m));
					}
				}
			}

			fprintf(fp, "\n");
		}
	}

	printf("Acceptance rate %f\n", (double)acc / (double)sample / (double)Nt);

	fclose(fp);

	free(lattice);
	free(nnp);
	free(nnm);

	double end = omp_get_wtime();
    double time_spent = (end - begin);
    printf("\nHa impiegato %.3f secondi\n", time_spent);

	return EXIT_SUCCESS;
}

int main(void)
{
	const unsigned long int seed1 = (unsigned long int)time(NULL);
	const unsigned long int seed2 = seed1 + 127;

	myrand_init(seed1, seed2);

	montecarlo(10, 5, 1e6);

	return EXIT_SUCCESS;
}
