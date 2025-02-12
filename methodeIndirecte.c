#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "methodeDirecte.h"
#include "methodeIndirecte.h"
#include "fonctionnalite.h"

#define N 3  // Taille de la matrice (N x N)
#define TOLERANCE 1e-6  // Tolérance pour la convergence
#define MAX_ITERATIONS 1000  // Nombre maximal d'itérations


// M�thode de Jacobi
void jacobi(double A[N][N], double B[N], double X[N]) {
    double X_new[N];
    int iterations = 0;
    double error;

    do {
        for (int i = 0; i < N; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * X[j];
            }
            for (int j = i + 1; j < N; j++) {
                sum2 += A[i][j] * X[j];
            }
            X_new[i] = (B[i] - sum1 - sum2) / A[i][i];
        }

        // Calculer l'erreur
        error = 0.0;
        for (int i = 0; i < N; i++) {
            error += fabs(X_new[i] - X[i]);
            X[i] = X_new[i];
        }

        iterations++;
    } while (error > TOLERANCE && iterations < MAX_ITERATIONS);

    if (iterations == MAX_ITERATIONS) {
        printf("La méthode de Jacobi n'a pas convergé après %d itérations.\n", MAX_ITERATIONS);
    }
}


 // M�thode de Gauss-Seidel
 void gaussSeidel(double A[N][N], double B[N], double X[N]) {
    int iterations = 0;
    double error;

    do {
        error = 0.0;
        for (int i = 0; i < N; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * X[j];
            }
            for (int j = i + 1; j < N; j++) {
                sum2 += A[i][j] * X[j];
            }
            double X_new = (B[i] - sum1 - sum2) / A[i][i];
            error += fabs(X_new - X[i]);
            X[i] = X_new;
        }

        iterations++;
    } while (error > TOLERANCE && iterations < MAX_ITERATIONS);

    if (iterations == MAX_ITERATIONS) {
        printf("La méthode de Gauss-Seidel n'a pas convergé après %d itérations.\n", MAX_ITERATIONS);
    }
}