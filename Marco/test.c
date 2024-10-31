#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<time.h>
#include"./include/random.h"

int main(void){
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;
    myrand_init(seed1, seed2);

    //scegli i valori di beta e L (ed eventualmente degli altri parametri)
    double a = myrand();    
    return EXIT_SUCCESS;
}
