#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/random.h"

#define SIZE 10
#define STRING_LENGTH 50

// magnetizzazione
double magn(int lattice[SIZE][SIZE][SIZE]) {
    long int rx, ry, rz, sum;

    sum = 0;
    for (rx = 0; rx < SIZE; rx++) {
        for (ry = 0; ry < SIZE; ry++) {
            for (rz = 0; rz < SIZE; rz++) {
                sum += lattice[rx][ry][rz];
            }
        }
    }

    return (double)sum / (double)(SIZE * SIZE * SIZE);
}

// energia
double energy(int lattice[SIZE][SIZE][SIZE]) {
    long int rx, ry, rz, tmp, sum;

    sum = 0;
    for (rx = 0; rx < SIZE; rx++) {
        for (ry = 0; ry < SIZE; ry++) {
            for (rz = 0; rz < SIZE; rz++) {
                tmp = (rx + 1) % SIZE;
                sum += -lattice[rx][ry][rz] * lattice[tmp][ry][rz];

                tmp = (ry + 1) % SIZE;
                sum += -lattice[rx][ry][rz] * lattice[rx][tmp][rz];

                tmp = (rz + 1) % SIZE;
                sum += -lattice[rx][ry][rz] * lattice[rx][ry][tmp];
            }
        }
    }

    return (double)sum / (double)(SIZE * SIZE * SIZE);
}

// metropolis 
int metropolis(int lattice[SIZE][SIZE][SIZE], 
               long int rx, long int ry, long int rz,
               double const * const restrict acc_prob) {
    long int tmp;
    int sumnn, acc = 0;

    // calcola solo il delta dell'energia
    sumnn = 0;

    // Neighbors in 3D
    tmp = (rx + 1) % SIZE;
    sumnn += lattice[tmp][ry][rz];

    tmp = (rx - 1 + SIZE) % SIZE;
    sumnn += lattice[tmp][ry][rz];

    tmp = (ry + 1) % SIZE;
    sumnn += lattice[rx][tmp][rz];

    tmp = (ry - 1 + SIZE) % SIZE;
    sumnn += lattice[rx][tmp][rz];

    tmp = (rz + 1) % SIZE;
    sumnn += lattice[rx][ry][tmp];

    tmp = (rz - 1 + SIZE) % SIZE;
    sumnn += lattice[rx][ry][tmp];

    sumnn *= lattice[rx][ry][rz];

    if (sumnn < 0) {
        lattice[rx][ry][rz] = -lattice[rx][ry][rz];
        acc = 1;
    } else {
        if (myrand() < acc_prob[sumnn]) {
            lattice[rx][ry][rz] = -lattice[rx][ry][rz];
            acc = 1;
        }
    }

    return acc;
}

// main
int main(void) {
    int i, lattice[SIZE][SIZE][SIZE];
    long int rx, ry, rz, rxaux, ryaux, rzaux, sample, iter, acc; 
    double beta, locE, locM;
    double acc_prob[7];

    char datafile[STRING_LENGTH];
    FILE *fp;

    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;

    beta = 0.2;
    sample = 1e7;
    strcpy(datafile, "./ising3d.dat");

    // inizializzazione generatore
    myrand_init(seed1, seed2);

    // inizializzazione reticolo
    for (rx = 0; rx < SIZE; rx++) {
        for (ry = 0; ry < SIZE; ry++) {
            for (rz = 0; rz < SIZE; rz++) {
                lattice[rx][ry][rz] = 1;
            }
        }
    }

    fp = fopen(datafile, "w");
    //fprintf(fp, "%d\n", SIZE);

    // si salva tutti i possibili delta dell'energia all'inizio
    for (i = 0; i < 7; i++) {
        acc_prob[i] = exp(-2.0 * beta * ((double)i));
    }

    // variabile che segna se è stato accettato il flip
    // ogni iterazione lui fa SIZE*SIZE*SIZE aggiornamenti
    acc = 0;
    for (iter = 0; iter < sample; iter++) {
        for (rx = 0; rx < SIZE; rx++) {
            for (ry = 0; ry < SIZE; ry++) {
                for (rz = 0; rz < SIZE; rz++) {
                    rxaux = (long int)((double)SIZE * myrand());
                    ryaux = (long int)((double)SIZE * myrand());
                    rzaux = (long int)((double)SIZE * myrand());

                    acc += metropolis(lattice, rxaux, ryaux, rzaux, acc_prob);
                }
            }
        }

        locE = energy(lattice);
        locM = magn(lattice);

        fprintf(fp, "%.12f %.12f\n", locE, locM);
    }

    fclose(fp);
    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double)(SIZE * SIZE * SIZE));

    return EXIT_SUCCESS;
}
