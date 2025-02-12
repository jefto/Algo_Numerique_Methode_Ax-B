#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "methodeDirecte.h"
#include "methodeIndirecte.h"
#include "fonctionnalite.h"

#define N 3  // Taille de la matrice (N x N)




// Résolution du système LY = B (substitution avant)
void substitutionAvant(double L[N][N], double B[N], double Y[N]) {
    for (int i = 0; i < N; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
        Y[i] /= L[i][i];
    }
}

// Résolution du système UX = Y (substitution arrière)
void substitutionArriere(double U[N][N], double Y[N], double X[N]) {
    for (int i = N - 1; i >= 0; i--) {
        X[i] = Y[i];
        for (int j = i + 1; j < N; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
}


// Fonction pour échanger deux lignes d'une matrice
void echangerLignes(double mat[N][N+1], int ligne1, int ligne2) {
    for (int j = 0; j <= N; j++) {
        double temp = mat[ligne1][j];
        mat[ligne1][j] = mat[ligne2][j];
        mat[ligne2][j] = temp;
    }
}

//Methode de gauss avec pivot
void gaussAvecPivot(double A[N][N+1]) {
    int i, j, k, maxIndex;
    double max, temp, ratio;

    for (i = 0; i < N; i++) {
        // Trouver le pivot
        max = fabs(A[i][i]);
        maxIndex = i;
        for (k = i + 1; k < N; k++) {
            if (fabs(A[k][i]) > max) {
                max = fabs(A[k][i]);
                maxIndex = k;
            }
        }

        // Échanger les lignes
        if (maxIndex != i) {
            for (j = 0; j <= N; j++) {
                temp = A[i][j];
                A[i][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }
        }

        // Élimination
        for (k = i + 1; k < N; k++) {
            ratio = A[k][i] / A[i][i];
            for (j = i; j <= N; j++) {
                A[k][j] -= ratio * A[i][j];
            }
        }
    }

    // Résolution du système
    double X[N];
    for (i = N - 1; i >= 0; i--) {
        X[i] = A[i][N];
        for (j = i + 1; j < N; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    // Afficher la solution
    printf("\nSolution du système :\n");
    for (i = 0; i < N; i++) {
        printf("X[%d] = %.2f\n", i, X[i]);
    }
}

void gaussSansPivot(double A[N][N+1]) {
    int i, j, k;
    double ratio;

    for (i = 0; i < N - 1; i++) {
        for (k = i + 1; k < N; k++) {
            ratio = A[k][i] / A[i][i];
            for (j = i; j <= N; j++) {
                A[k][j] -= ratio * A[i][j];
            }
        }
    }

    // Résolution du système
    double X[N];
    for (i = N - 1; i >= 0; i--) {
        X[i] = A[i][N];
        for (j = i + 1; j < N; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    // Afficher la solution
    printf("\nSolution du système :\n");
    for (i = 0; i < N; i++) {
        printf("X[%d] = %.2f\n", i, X[i]);
    }
}

void gaussJordan(double A[N][N], double B[N], double X[N]) {
    double matriceAugmentee[N][N+1];

    // Création de la matrice augmentée [A|B]
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matriceAugmentee[i][j] = A[i][j];
        }
        matriceAugmentee[i][N] = B[i];
    }

    printf("Matrice augmentée avant élimination de Gauss-Jordan :\n");
    afficherMatriceAugmentee(matriceAugmentee);

    // Élimination de Gauss-Jordan
    for (int k = 0; k < N; k++) {
        // Normalisation de la ligne k
        double pivot = matriceAugmentee[k][k];
        for (int j = 0; j <= N; j++) {
            matriceAugmentee[k][j] /= pivot;
        }

        // Élimination des autres lignes
        for (int i = 0; i < N; i++) {
            if (i != k) {
                double facteur = matriceAugmentee[i][k];
                for (int j = 0; j <= N; j++) {
                    matriceAugmentee[i][j] -= facteur * matriceAugmentee[k][j];
                }
            }
        }
    }

    printf("\nMatrice augmentée après élimination de Gauss-Jordan :\n");
    afficherMatriceAugmentee(matriceAugmentee);

    // Extraction des solutions
    for (int i = 0; i < N; i++) {
        X[i] = matriceAugmentee[i][N];
    }

    // Afficher la solution
    printf("\nSolution du système :\n");
    afficherVecteur(X);
}


//Methode de thomas
void thomas(double a[], double b[], double c[], double d[], double x[], int n) {
    double c_prime[n];
    double d_prime[n];

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double m = 1.0 / (b[i] - a[i-1] * c_prime[i-1]);
        c_prime[i] = c[i] * m;
        d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) * m;
    }

    // Back substitution
    x[n-1] = d_prime[n-1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
}

// Décomposition LU Doolittle
void doolittle(double A[N][N], double L[N][N], double U[N][N]) {
    // Initialisation de L et U
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                L[i][j] = 1.0;  // Diagonale de L remplie de 1
            } else {
                L[i][j] = 0.0;
            }
            U[i][j] = 0.0;
        }
    }

    // Décomposition
    for (int i = 0; i < N; i++) {
        // Calcul de U
        for (int k = i; k < N; k++) {
            double somme = 0.0;
            for (int j = 0; j < i; j++) {
                somme += L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - somme;
        }

        // Calcul de L
        for (int k = i + 1; k < N; k++) {
            double somme = 0.0;
            for (int j = 0; j < i; j++) {
                somme += L[k][j] * U[j][i];
            }
            L[k][i] = (A[k][i] - somme) / U[i][i];
        }
    }
}

void luCrout(double A[N][N], double L[N][N], double U[N][N]) {
    int i, j, k;

    // Initialiser L et U
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i >= j) {
                U[i][j] = 0;
            } else {
                L[i][j] = 0;
            }
        }
    }

    // Décomposition LU Crout
    for (j = 0; j < N; j++) {
        for (i = 0; i <= j; i++) {
            L[i][j] = A[i][j];
            for (k = 0; k < i; k++) {
                L[i][j] -= L[i][k] * U[k][j];
            }
        }
        for (i = j; i < N; i++) {
            U[j][i] = A[j][i] / L[j][j];
            for (k = 0; k < j; k++) {
                U[j][i] -= (L[j][k] * U[k][i]) / L[j][j];
            }
        }
    }
}
