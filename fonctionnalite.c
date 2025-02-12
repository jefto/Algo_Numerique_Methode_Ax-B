#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "methodeDirecte.h"
#include "methodeIndirecte.h"
#include "fonctionnalite.h"

#define N 3  // Taille de la matrice (N x N)
#define TOLERANCE 1e-6  // Tolérance pour la convergence
#define MAX_ITERATIONS 1000  // Nombre maximal d'itérations

double A[N][N+1] = {
    {-3, -1, 2, -11},
    {2, 1, -1, 8},
    {-2, 1, 2, -3}
    
};

double As[N][N] = {
    {-3, -1, 2},
    {2, 1, -1},
    {-2, 1, 2}
    
};

//methode de thomas
double A_tridiag[N][N] = {
    {2, -1, 0},
    {-1, 2, -1},
    {0, -1, 2}
};
double B_tridiag[N] = {1, 2, 3};

double B[N] = {-11, 8, -3};

// Fonction pour afficher une matrice
void afficherMatrice(double matrice[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.2f\t", matrice[i][j]);
        }
        printf("\n");
    }
}

void afficherMatriceAugmentee(double matrice[N][N+1]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= N; j++) {
            printf("%.2f\t", matrice[i][j]);
        }
        printf("\n");
    }
}

// Fonction pour afficher un vecteur
void afficherVecteur(double vecteur[N]) {
    for (int i = 0; i < N; i++) {
        printf("X[%d] = %.6f\n", i, vecteur[i]);
    }
}

void menuRelancer(){
    printf("\nVoulez vous utiliser une autre méthode?\n");
    printf("1. Oui\n");
    printf("2. Non\n");

    printf("\nVeuillez entrer votre choix (1-2) : ");

    int choix;
    scanf("%d", &choix);
    switch (choix)
    {
    case 1:
        afficherMenuPrincipal();
        break;
    case 2:
        printf("Au revoir !\n");
        break;
    default:
        printf("Choix entrer invalide. Veuillez entrer un nombre entier entre 1 et 2.\n");
        menuRelancer();
        break;
    }

}

void afficherMenuPrincipal(){
    printf("--------------------------\n");
    printf("--------------------------\n");
    printf("|||||||||| MENU ||||||||||\n");
    printf("--------------------------\n");
    printf("--------------------------\n");

    printf("\n 1. Méthode Directe\n");
    printf(" 2. Méthode Indirecte\n");
    printf("\n 0. Quitter\n");

    printf("\nVeuillez entrer votre choix (0-2) : ");

    choixMenuPrincipal();
}

void afficherMenuMethodeDirecte(){
    printf("|||||||||| METHODES DIRECTS ||||||||||\n");
    printf("\n 1. Gauss avec Pivot\n");
    printf(" 2. Gauss sans Pivot\n");
    printf(" 3. Gauss-Jordan\n");
    printf(" 4. Thomas\n");
    printf(" 5. Doolittle\n");
    printf(" 6. Lu Crout\n");
    printf("\n 0. Retour\n");

    printf("\nVeuillez entrer votre choix (0-5) : ");

    choixMenuMethodeDirecte();
}

void afficherMethodeMenuIndirecte(){
    printf("|||||||||| METHODES INDIRECTS ||||||||||\n");
    printf("\n 1. Jacobi\n");
    printf(" 2. Gauss-Seidel\n");
    printf("\n 0. Retour\n");

    printf("\nVeuillez entrer votre choix (0-2) : ");

    choixMenuMethodeIndirecte();
}

void choixMenuPrincipal(){
    int choix;

    scanf("%d", &choix);
    switch (choix) {
        case 1:
            afficherMenuMethodeDirecte();
            break;
        case 2:
            afficherMethodeMenuIndirecte();
            break;
        case 0:
            printf("Au revoir !\n");
            break;
        default:
            printf("Choix entrer invalide. Veuillez entrer un nombre entier entre 0 et 2.\n");
            choixMenuPrincipal();
            break;
    }
}

void choixMenuMethodeDirecte(){
    int choix;
    double X[N] = {0, 0, 0};  // Vecteur solution
    double L[N][N], U[N][N];  // Matrices L et U
    double Y[N];  // Vecteur Y

    // Coefficients pour la méthode de Thomas
    /*
    double a[N-1] = {-3, -2};  // Coefficients sous-diagonaux
    double b[N] = {2, -1, 2};  // Coefficients diagonaux
    double c[N-1] = {1, 1};  // Coefficients super-diagonaux
    double d[N] = {8, -11, -3};  // Termes constants
    */
    

    double a[N-1] = {-1, -1};  // Coefficients sous-diagonaux
    double b[N] = {2, 2, 2};  // Coefficients diagonaux
    double c[N-1] = {-1, -1};  // Coefficients super-diagonaux
    double d[N] = {1, 2, 3};  // Termes constants

    // Extraire la matrice A et le vecteur B de la matrice augmentée(Methode de gaussJordan)
    double A_simple[N][N];
    double B_simple[N];


    scanf("%d", &choix);
    switch (choix) {
        case 1:
            printf("Matrice initiale :\n");
            afficherMatriceAugmentee(A);

            gaussAvecPivot(A);
            menuRelancer();
            break;
        case 2:
            printf("Matrice initiale :\n");
            afficherMatriceAugmentee(A);

            gaussSansPivot(A);
            menuRelancer();
            break;
        case 3:
             
             for (int i = 0; i < N; i++) {
                 for (int j = 0; j < N; j++) {
                     A_simple[i][j] = A[i][j];
                 }
                 B_simple[i] = A[i][N];
             }
 
             gaussJordan(A_simple, B_simple, X);
 
             printf("\nSolution du système :\n");
             for (int i = 0; i < N; i++) {
                 printf("X[%d] = %.2f\n", i, X[i]);
             }
             menuRelancer();
            break;
        case 4:
           thomas(a, b, c, d, X, N);
            printf("\nSolution du système :\n");
            for (int i = 0; i < N; i++) {
                printf("X[%d] = %.2f\n", i, X[i]);
            }
            menuRelancer();
            break;
        case 5:
            // Décomposition LU doolittle
            doolittle(As, L, U);

            printf("Matrice L :\n");
            afficherMatrice(L);
            printf("\nMatrice U :\n");
            afficherMatrice(U);

            // Résolution du système LY = B
            substitutionAvant(L, B, Y);
            printf("\nVecteur Y :\n");
            afficherVecteur(Y);

            // Résolution du système UX = Y
            substitutionArriere(U, Y, X);
            printf("\nSolution X :\n");
            afficherVecteur(X);
            menuRelancer();
            break;
        case 6:
            luCrout(As, L, U);

            printf("Matrice L :\n");
            afficherMatrice(L);
            printf("\nMatrice U :\n");
            afficherMatrice(U);

            // Résolution du système LY = B
            substitutionAvant(L, B, Y);
            printf("\nVecteur Y :\n");
            afficherVecteur(Y);

            // Résolution du système UX = Y
            substitutionArriere(U, Y, X);
            printf("\nSolution X :\n");
            afficherVecteur(X);
            menuRelancer();
        break;
            break;
        case 0:
            afficherMenuPrincipal();
            break;
        default:
            printf("Choix entrer invalide. Veuillez entrer un nombre entier entre 0 et 5.\n");
            afficherMenuMethodeDirecte();
            break;
    }
}

void choixMenuMethodeIndirecte(){
    int choix;
    double X[N] = {0, 0, 0};    // Solution initiale

    scanf("%d", &choix);
    switch (choix) {
        case 1:
            // Résolution par Jacobi
            jacobi(As, B, X);

            // Affichage de la solution
            printf("\nSolution du système :\n");
            afficherVecteur(X);
            menuRelancer();
            break;
        case 2:
            // Résolution par Gauss-Seidel
            gaussSeidel(As, B, X);

            // Affichage de la solution
            printf("\nSolution du système :\n");
            afficherVecteur(X);
            menuRelancer();
            break;
        case 0:
            afficherMenuPrincipal();
            break;
        default:
            printf("Choix entrer invalide. Veuillez entrer un nombre entier entre 0 et 2.\n");
            afficherMethodeMenuIndirecte();
            break;
    }
}
