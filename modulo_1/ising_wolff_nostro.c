#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include <omp.h>
#include"./include/geometry.h"
#include"./include/random.h"
#include"./include/get_array.h"

#define DIM 2 // dimensionality
#define STRING_LENGTH 50

// magnetization per site
double magn(int const * const restrict lattice, long int volume)
  {
  long int r, sum;

  sum=0;
  for(r=0; r<volume; r++)
     {
     sum+=lattice[r];
     }

  return (double) sum / (double) volume;
  }


void init_neighbours_mio(long int *nn, int L, int q, void (*nearest)(int, int, int *, int *, int)){
    int i, j, nnx[q], nny[q];

    for (i = 0; i < L*L; i++){
        nearest(i/L, i%L, nnx, nny, L);
        for (j = 0; j < q; j++){
            nn[i * q + j] = nnx[j] * L + nny[j];
        }
    }

}


// non-recursive construction of the culter
void build_cluster_norec(int const * const restrict lattice, 
                         long int r, 
                         int * restrict occup, 
                         long int * restrict pointtoocc, 
                         long int * restrict clustersize,
                         long int const * const restrict nn,
                         long int volume,
                         double prob,
                         int q)
  {
  #define nn_pos(rr, ii) ((rr) * q + (ii))
  (void) r; // just to avoid warnings
  int i;
  long int index, r1;
  long int oldcs, oldcsnew;

  oldcs=0; // starting value, with *clustersize=1

  // if first neighbors have the same orientation and are not occupied
  // they are added to the cluster with probability prob

  while(*clustersize>oldcs) // this means that there are sites recently added, whose neighbors has not been checked yet, so we check them
       { 
       oldcsnew=*clustersize;

       for(index=oldcs; index<oldcsnew; index++)
          {
          r1=pointtoocc[index];

          for(i=0; i<q; i++)
             {
             if(occup[nn[nn_pos(r1, i)]]==0 && lattice[r1]*lattice[nn[nn_pos(r1, i)]]==1)
               {
               if(myrand()<prob)
                 {
                 occup[nn[nn_pos(r1, i)]]=1;
                 pointtoocc[*clustersize]=nn[nn_pos(r1, i)];
                 (*clustersize)++;
                 }
               }
             } 
          }

       oldcs=oldcsnew;
       }
  }


// main
int main()
    {
    
    char *modello_values[] = {"ising2d_sq_cluster"}; //{"ising2d_tri_cluster", "ising2d_sq_cluster", "ising2d_hex_cluster"};
    int num_modelli = sizeof(modello_values) / sizeof(modello_values[0]);
    
    int num_L = 1;
    int L_array[] = {100};

    int num_beta = 40;
    
    for (int mm = 0; mm < num_modelli; mm++) {       
        #pragma omp parallel for collapse(2) shared(L_array, modello_values, num_beta) schedule(dynamic, 1)  // collapse the loops and define private variables
        for (int ii = 0; ii < num_beta; ii++){
            for (int jj = 0; jj < num_L; jj++){
              int L = L_array[jj];
              int i, *lattice, *occup;
              long int r, volume, sample, iter, clustersize; 
              long int *nn, *pointtoocc;
              double locE, locM, prob;

              int q;
              void (*nearest)(int, int, int *, int *, int);
              double (*energy)(int *restrict, int, int, double);
              
              char *modello = modello_values[mm];
              choose_geometry(modello, &nearest, &energy, &q);
              int nnx[q], nny[q];
    
              double beta = assign_beta_close(modello, ii, num_beta, L_array[jj]);
              
            
              sample = 1e6;
              char datafile[STRING_LENGTH];
              sprintf(datafile, "./%s/L%d_beta%.5f.dat", modello, L, beta); // file name initialized with a string
    
              FILE *fp;

              const unsigned long int seed1=(unsigned long int) time(NULL) + omp_get_thread_num();
              const unsigned long int seed2=seed1+127;

              // initialize random number generator
              myrand_init(seed1, seed2);

              // compute the volume
              volume=1;
              for(i=0; i<DIM; i++)
              {
              volume*=L;
              }

              // allocate the lattice (lexicographic order)
              // and next neighbors: nnp[dirgeo(r, i, volume)]= next neighbor in positive "i" direction of site r 
              lattice=(int *)malloc((unsigned long int)(volume)*sizeof(int));
              if(lattice == NULL)
              {
              fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
              //return EXIT_FAILURE;
              }
              nn=(long int *)malloc((unsigned long int)(q*volume)*sizeof(long int));
              if(nn == NULL){
              fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
              //return EXIT_FAILURE;
              }

              // this structure will be used to keep trak of the occupied sites 
              // while building the cluster
              // 0 = free, 1=occupied
              occup=(int *)malloc((unsigned long int)(volume)*sizeof(int));
              if(occup== NULL)
              {
              fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
              //return EXIT_FAILURE;
              }

              // the first "clustesize" entries will point to the cluster sites
              pointtoocc=(long int *)malloc((unsigned long int)(volume)*sizeof(long int));
              if(pointtoocc== NULL)
              {
              fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
              //return EXIT_FAILURE;
              }

              // initialize nn
              init_neighbours_mio(nn, L, q, nearest);

              // initialize lattice to ordered start
              for(r=0; r<volume; r++)
              {
              lattice[r]=1;
              }

              // open data file
              fp=fopen(datafile, "w");
              if(fp==NULL)
              {
              fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
              //return EXIT_FAILURE;
              }

              fprintf(fp, "m, E, beta = %.10f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %ld\n", beta, L, 0, 1, sample); 
              // probability of addition to the cluster
              prob=1.0-exp(-2.0*beta);

              for(iter=0; iter<sample; iter++)
              {
                for(r=0; r<volume; r++) 
                    {
                    occup[r]=0;
                    }
                clustersize=0;

                r=(int)((double)volume*myrand());
                occup[r]=1; // r is set as occupied
                pointtoocc[clustersize]=r; // a pointer to "r" is added in position "clustersize"
                clustersize++;

                //build_cluster_rec(lattice, r, occup, pointtoocc, &clustersize, nnp, nnm, volume, prob);
                build_cluster_norec(lattice, r, occup, pointtoocc, &clustersize, nn, volume, prob, q);

                // flip the cluster
                for(r=0; r<clustersize; r++)
                    {
                    lattice[pointtoocc[r]]=-lattice[pointtoocc[r]];
                    }

                locE=energy(lattice, volume, L, beta);
                locM=magn(lattice, volume);

                fprintf(fp, "%.10f, %.10f\n", locM, locE);
              }

              // close datafile
              fclose(fp);

              free(lattice);
              free(occup);
              free(pointtoocc);
              free(nn);
          
          }
        }
    }
    
    return EXIT_SUCCESS;
}
