#ifndef FONCTIONNALITE_H
#define FONCTIONNALITE_H

#define N 3

    extern double A[N][N+1];

    extern double As[N][N];

    extern double B[N];

    extern double A_tridiag[N][N];

    double B_tridiag[N];

    void afficherVecteur(double vecteur[N]);
    void afficherMatrice(double matrice[N][N]);
    void afficherMatriceAugmentee(double matrice[N][N+1]);

    void afficherMenuPrincipal();
    void afficherMenuMethodeDirecte();
    void afficherMethodeMenuIndirecte();
    void choixMenuPrincipal();
    void choixMenuMethodeDirecte();
    void choixMenuMethodeIndirecte();
#endif
