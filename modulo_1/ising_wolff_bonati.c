//OPEN ################################################################
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
//CLOSE ################################################################

// energy per site
double energy(int const * const restrict lattice, 
              long int const * const restrict nnp, 
              long int volume)
  {
  long int r, sum;
  int i;

  sum=0;
  for(r=0; r<volume; r++)
     {
     for(i=0; i<DIM; i++)
        {
        sum+=-lattice[r]*lattice[nnp[dirgeo(r, i, volume)]];
        }
     }

  return (double) sum / (double) volume;
  }

//OPEN ################################################################
// non-recursive construction of the culter
void build_cluster_norec(int const * const restrict lattice, 
                         long int r, 
                         int * restrict occup, 
                         long int * restrict pointtoocc, 
                         long int * restrict clustersize,
                         long int const * const restrict nnp, //DIVERSO
                         long int const * const restrict nnm, // DIVERSO
                         long int volume,
                         double prob)
  {
  (void) r; // just to avoid warnings
  int i;
  long int index, r1;
  long int oldcs, oldcsnew;

  oldcs=0; // starting value, with *clustersize=1

  // if first neighbors have the same orientation and are not occupied
  // they are added to the cluster with probability prob
//CLOSE ################################################################
  while(*clustersize>oldcs) // this means that there are sites recently added, whose neighbors has not been checked yet, so we check them
       { 
       oldcsnew=*clustersize;

       for(index=oldcs; index<oldcsnew; index++)
          {
          r1=pointtoocc[index];

          for(i=0; i<DIM; i++)
             {
             // forward
             if(occup[nnp[dirgeo(r1, i, volume)]]==0 && lattice[r1]*lattice[nnp[dirgeo(r1, i, volume)]]==1)
               {
               if(myrand()<prob)
                 {
                 occup[nnp[dirgeo(r1, i, volume)]]=1;
                 pointtoocc[*clustersize]=nnp[dirgeo(r1, i, volume)];
                 (*clustersize)++;
                 }
               }
        
             // backward
             if(occup[nnm[dirgeo(r1, i, volume)]]==0 && lattice[r1]*lattice[nnm[dirgeo(r1, i, volume)]]==1) 
               {
               if(myrand()<prob)
                 {
                 occup[nnm[dirgeo(r1, i, volume)]]=1;
                 pointtoocc[*clustersize]=nnm[dirgeo(r1, i, volume)];
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
    
    int num_beta = 40;
    
    #pragma omp parallel for shared(num_beta) schedule(dynamic, 1) // collapse the loops and define private variables
    for (int i = 0; i < num_beta; i++){
        int L = 100;
        double beta = assign_beta_close("ising2d_sq", i, num_beta, L);
        int i, *lattice, *occup;
        long int r, volume, sample, iter, clustersize; 
        long int *nnp, *nnm, *pointtoocc;
        double locE, locM, prob;
        
        sample = 1e6;
        char datafile[STRING_LENGTH];
        sprintf(datafile, "./%s/L%d_beta%.5f.dat", "bonati", L, beta);
        FILE *fp;

        const unsigned long int seed1=(unsigned long int) time(NULL);
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
        nnp=(long int *)malloc((unsigned long int)(DIM*volume)*sizeof(long int));
        if(nnp == NULL){
        fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        //return EXIT_FAILURE;
        }
        nnm=(long int *)malloc((unsigned long int)(DIM*volume)*sizeof(long int));
        if(nnm == NULL){
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

        // initialize nnp and nnm
        init_neighbors(nnp, nnm, L, DIM);

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

        fprintf(fp, "m, E, beta = %f, L = %d, iterations = %d, iter_bet_meas = %d, num_measures = %ld\n", beta, L, 0, 1, sample); 
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
            build_cluster_norec(lattice, r, occup, pointtoocc, &clustersize, nnp, nnm, volume, prob);

            // flip the cluster
            for(r=0; r<clustersize; r++)
                {
                lattice[pointtoocc[r]]=-lattice[pointtoocc[r]];
                }

            locE=energy(lattice, nnp, volume);
            locM=magn(lattice, volume);

            fprintf(fp, "%.12f %.12f\n", locM, locE);
        }

        // close datafile
        fclose(fp);

        free(lattice);
        free(occup);
        free(pointtoocc);
        free(nnp);
        free(nnm);
    
    }
    
    return EXIT_SUCCESS;
}
