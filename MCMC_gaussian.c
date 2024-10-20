#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_ITER 10000  // Numero di iterazioni

// Funzione per generare un numero casuale tra un range
double rand_uniform(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

int main() {
    double x0 = 5.0;
    double delta = 0.1;
    double xk = x0;
    double xkp1;
    double y, r;
    double history[N_ITER];
    
    // Inizializzazione del generatore di numeri casuali
    srand(time(NULL));
    
    // Simulazione del loop
    for (int k = 0; k < N_ITER; k++) {
        // Seleziona x nell'intervallo (xk - delta, xk + delta)
        double x = rand_uniform(xk - delta, xk + delta);
        
        // Calcola y
        y = -0.5 * x * x + 0.5 * xk * xk;
        
        // Condizione per aggiornare xkp1
        if (y >= 0) {
            xkp1 = x;
        } else {
            // Seleziona r ∈ [0, 1) con distribuzione uniforme
            r = rand_uniform(0, 1);
            if (r <= fmin(1, exp(y))) {
                xkp1 = x;
            } else {
                xkp1 = xk;
            }
        }
        
        // Salva il nuovo valore di x
        history[k] = xkp1;
        
        // Aggiorna xk
        xk = xkp1;
    }
    
    // Stampa i risultati in un file per la grafica
    FILE *file = fopen("output.txt", "w");
    if (file == NULL) {
        printf("Errore nell'apertura del file!\n");
        return 1;
    }
    
    for (int i = 0; i < N_ITER; i++) {
        fprintf(file, "%d %f\n", i, history[i]);
    }
    
    fclose(file);
    
    // Istruzioni per graficare con gnuplot
    printf("Per visualizzare il risultato, usa il seguente comando in gnuplot:\n");
    printf("plot 'output.txt' with lines title 'x(k) in funzione di k'\n");

    return 0;
}

