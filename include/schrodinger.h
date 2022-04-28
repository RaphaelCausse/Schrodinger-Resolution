#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include <common.h>

#define NB_CMDS     6      //Nombre de commandes gnuplot
#define NB_PTS      100    //Nombre de subdivision de [xmin, xmax]
#define HBARC       197    //hbar*c en eV*nm (voir main.c pour les explications)

typedef struct _psi {
    double psi;
    double dpsi;
} _Psi;

typedef struct _solve {
    double *dpsi;
    double *E;
} _Solve;


/********** Utils **********/

/*  
    Affiche un tableau de taille n.
*/
void display_array(double *array, const size_t n, const char *name);


/********** Resolution de l'equation de Schrodinger independante du temps **********/

/*
*/
void store_data(double xdata[], double psidata[], double x, double y[], double f[], void *params_ptr);

/*
*/
void write_data(FILE *data, double xdata[], double psidata[]);

/*
    Plot graphiques de donnees avec gnuplot.
    mode 0: puit de potentiel infini
*/
void plot_gnuplot(int mode);

/*
    Methode d'Euler pour resoudre une eq diff.
*/
void euler_method(double x, double y[], double f[], void *params_ptr);

/*
    Resolution du systeme avec condition de normalisation, en utilisant la methode d'Euler.
*/
void solve_euler(_Psi *p, double psi, double dpsi, double m, double E, double V, double xmin, double xmax);

/*
    Resolution dans le cas du puit de potentiel infini.
*/
void puit_potentiel_infini(double m, double L, double E);

#endif