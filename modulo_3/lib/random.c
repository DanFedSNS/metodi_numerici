#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "../include/pcg32min.h"
#include "../include/random.h"

// Thread-local RNG state
#pragma omp threadprivate(pcg32_random_state) 
pcg32_random_t pcg32_random_state;

// Initialize RNG for a thread
void myrand_init(unsigned long int initstate, unsigned long int initseq) {
    pcg32_srandom_r(&pcg32_random_state, (uint64_t) initstate, (uint64_t) initseq);
}

// Generate a random number in [0,1)
double myrand(void) {
    return (double) pcg32_random_r(&pcg32_random_state) / (pow(2.0, 32.0));
}
