#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <stdbool.h>
#include<time.h>

int L = 100; // Replace with your desired dimension size

void first(){
    double a = cos(1.23);
}

void second(){
    double a = 2;
    double b = a;
}

void time_func(void (*func)(void), int rep, int caso){
    double time_spent;

    clock_t begin = clock();
    
    for (int j = 0; j < rep; j++){
        func();
    }

    clock_t end = clock();
    
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nCaso %d ha impiegato %.3f secondi\n", caso, time_spent);
}

int main(void){
    int ripetizioni = 1e9;
    
    time_func(first, ripetizioni, 1);
    time_func(second, ripetizioni, 2);

    clock_t begin = clock();
    
    clock_t end = clock();  
}