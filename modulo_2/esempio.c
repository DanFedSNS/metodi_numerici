#include <stdio.h>
#include <stdlib.h>
#include <arpack/arpack.h>

void matvec(const double *A, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i * n + j] * x[j];
        }
    }
}

void test_arpack() {
    int n = 100;  // Small matrix size for debugging
    int nev = 1; // Fewer eigenvalues
    int ncv = 5; // Reduced Lanczos vectors
    int ido = 0;
    char bmat = 'I';
    char which[2] = "SM";
    double tol = 1e-6;
    double *resid = malloc(n * sizeof(double));
    double *v = malloc(n * ncv * sizeof(double));
    int ldv = n;
    int iparam[11] = {0};
    int ipntr[14];
    double *workd = malloc(3 * n * sizeof(double));
    int lworkl = ncv * (ncv + 8); // Correct size
    double *workl = malloc(lworkl * sizeof(double));
    int info = 0;

    for (int i = 0; i < n; i++) {
        resid[i] = 0.1;
    }

    // Symmetric matrix
    double *A = malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = (i == j) ? 2.0 : (abs(i - j) == 1 ? 1.0 : 0.0);
        }
    }

    iparam[0] = 1;
    iparam[2] = 100; // Smaller max iterations
    iparam[6] = 1;

    while (ido != 99) {
        dsaupd_c(&ido, &bmat, n, which, nev, tol, resid, ncv, v, ldv,
                 iparam, ipntr, workd, workl, lworkl, &info);



        if (info != 0) {
            printf("Error during dsaupd_c: info = %d\n", info);
            break;
        }

        if (ido == 1) {
            const double *x = &workd[ipntr[0] - 1];
            double *y = &workd[ipntr[1] - 1];
            matvec(A, x, y, n);
        }
    }

    free(resid);
    free(v);
    free(workd);
    free(workl);
    free(A);

    if (info == 0) {
        printf("Eigenvalues converged successfully.\n");
    }
    for (int i = 0; i<nev; i++){
        printf("Eigenvalue %d: %f\n", i + 1, v[i * n + i]);
    }
}

int main() {
    test_arpack();
    return 0;
}
