#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include <common.h>

#define NB_CMDS     6      //Number of gnuplot commands
#define NB_PTS      1000   //Number of subdivision of interval [xmin, xmax]
#define HBARC       197    //hbar*c in eV*nm (see main.c for explanations)

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
    Display an array of length n.
*/
void display_array(double *array, const size_t n, const char *name);


/********** Solving Time-Independent Schrodinger equation **********/

/*
    
*/
void store_data(double xdata[], double psidata[], double x, double y[], double f[], void *params_ptr);

/*
*/
void write_data(FILE *data, double xdata[], double psidata[]);

/*
    Plot data graphs using gnuplot
    mode 0: infinte potential well
*/
void plot_gnuplot(int mode);

/*
    Euler's method to solve differential equation.
*/
void euler_method(double x, double y[], double f[], void *params_ptr);

/*
    Solve ODE system with normalisaiton condition, using Euler's method.
*/
void solve_euler(_Psi *p, double psi, double dpsi, double m, double E, double V, double xmin, double xmax);

/*
    Solve case of infinite potential well.
*/
void infinite_potential_well(double m, double L, double E);

#endif