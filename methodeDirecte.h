#ifndef METHODE_DIRECTE_H
#define METHODE_DIRECTE_H

#define N 3

    void echangerLignes(double mat[N][N+1], int ligne1, int ligne2);

    void gaussAvecPivot(double mat[N][N+1]);
    
    void gaussSansPivot(double mat[N][N+1]);

    void gaussJordan(double A[N][N], double B[N], double X[N]);

    void doolittle(double A[N][N], double L[N][N], double U[N][N]);

        void substitutionAvant(double L[N][N], double B[N], double Y[N]) ;

        void substitutionArriere(double U[N][N], double Y[N], double X[N]);

    void thomas(double a[], double b[], double c[], double d[], double x[], int n);

    void luCrout(double A[N][N], double L[N][N], double U[N][N]);
#endif
