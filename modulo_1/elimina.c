#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "./include/random.h"

#include <omp.h>


void nearest_sq(int rx, int ry, int *resx, int *resy, int L){  //nn reticolo quadrato
    resx[0] = (rx - 1 + L) % L;     //+ L serve a evitare negativi
    resx[1] = (rx + 1) % L;
    resx[2] = rx;
    resx[3] = rx;

    resy[0] = ry;
    resy[1] = ry;
    resy[2] = (ry - 1 + L) % L;
    resy[3] = (ry + 1) % L;
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



