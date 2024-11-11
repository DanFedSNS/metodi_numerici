#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"

#define SIZE 10
#define STRING_LENGTH 50

// magnetizzazione
double magn(int lattice[SIZE][SIZE])
  {
  long int rx, ry, sum;

  sum=0;
  for(rx=0; rx<SIZE; rx++)
     {
     for(ry=0; ry<SIZE; ry++)
        {
        sum+=lattice[rx][ry];
        }
     }

  return (double) sum / (double) (SIZE*SIZE);
  }


// energia
double energy(int lattice[SIZE][SIZE]) 
  {
  long int rx, ry, tmp, sum;

  sum=0;
  for(rx=0; rx<SIZE; rx++)
     {
     for(ry=0; ry<SIZE; ry++)
        {
        tmp=(rx+1)%SIZE;
        sum+=-lattice[rx][ry]*lattice[tmp][ry];
       
        tmp=(ry+1)%SIZE;
        sum+=-lattice[rx][ry]*lattice[rx][tmp];
        }
     }

  return (double) sum / (double) (SIZE*SIZE);
  }


// metropolis 
int metropolis(int lattice[SIZE][SIZE], 
               long int rx,
               long int ry,
               double const * const restrict acc_prob)
  {
  long int tmp;
  int sumnn, acc=0;

  // calcola solo il delta dell'energia
  sumnn=0;
  
  tmp=(rx+1)%SIZE;
  sumnn+=lattice[tmp][ry];

  tmp=(rx-1);
  if(tmp<0)
    {
    tmp=SIZE-1;
    }
  sumnn+=lattice[tmp][ry];

  tmp=(ry+1)%SIZE;
  sumnn+=lattice[rx][tmp];

  tmp=(ry-1);
  if(tmp<0)
    {
    tmp=SIZE-1;
    }
  sumnn+=lattice[rx][tmp];

  sumnn*=lattice[rx][ry];

  if(sumnn<0)
    {
    lattice[rx][ry]=-lattice[rx][ry];
    acc=1;
    }
  else
    {
    if(myrand()<acc_prob[sumnn])  // remember that acc_prob[i]=exp(-2.0*beta*((double)i));
      {
      lattice[rx][ry]=-lattice[rx][ry];
      acc=1;
      }
    }

  return acc;
  }


// main
int main(void)
    {
    int i, lattice[SIZE][SIZE];
    long int rx, ry, rxaux, ryaux, sample, iter, acc; 
    double beta, locE, locM;
    double acc_prob[5];
  
    char datafile[STRING_LENGTH];
    FILE *fp;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    beta=0.4;
    sample=1e6;
    strcpy(datafile, "./ising2d.dat");


    //inizializzazione generatore
    myrand_init(seed1, seed2);

    //inizializzazione reticolo
    for(rx=0; rx<SIZE; rx++)
       {
       for(ry=0; ry<SIZE; ry++)
          {
          lattice[rx][ry]=1;
          }
       }

    fp=fopen(datafile, "w");

    //si salva tutti i possibili delta dell'energia all'inizio
    for(i=0; i<5; i++)
       {
       acc_prob[i]=exp(-2.0*beta*((double)i));
       }
    
    //variabile che segna se è stato accettato il flip
    //ogni iterazione lui fa SIZE*SIZE aggiornamenti
    acc=0;
    for(iter=0; iter<sample; iter++)
       {
       for(rx=0; rx<SIZE; rx++)
          {
          for(ry=0; ry<SIZE; ry++)
             {
             rxaux=(long int)((double)SIZE * myrand());
             ryaux=(long int)((double)SIZE * myrand());

             acc+=metropolis(lattice, rxaux, ryaux, acc_prob);
             }
          }

       locE=energy(lattice);
       locM=magn(lattice);

       fprintf(fp, "%.12f %.12f\n", locE, locM);
       }

    fclose(fp);
    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) (SIZE*SIZE));

    return EXIT_SUCCESS;
    }


