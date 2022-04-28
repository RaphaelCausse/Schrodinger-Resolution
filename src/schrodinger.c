#include <schrodinger.h>

/********** Utils **********/

void display_array(double *array, const size_t n, const char *name) {
    printf("%s:\n\t[", name);
    for (int i = 0; i < n; i++) {
        printf(" %lf ", array[i]);
        if (i != n-1) printf(";");
    }
    printf("]\n");
}

/********** Resolution de l'equation de Schrodinger independante du temps **********/

void store_data(double xdata[], double psidata[], double x, double y[], double f[], void *params_ptr) {
    double *params = (double *)params_ptr;
    double xmax = params[0];
    double h = params[1];
    double k = params[2];
    double E = params[3];
    double V = params[4];

    int i = 0;
    while (x < xmax && i < NB_PTS) {
        f[1] = (-k)*(E-V)*y[0];
        f[0] = f[0] + h*f[1];
        y[1] = f[0];
        y[0] = y[0] + h*y[1];
        psidata[i] = y[0];
        xdata[i] = x;
        x += h;
        i++;
    }
}

void write_data(FILE *data, double xdata[], double psidata[]) {
    for (int i = 0; i < NB_PTS; i++) {
        fprintf(data, "%lf %lf\n", xdata[i], psidata[i]);
    }
    fprintf(data, "\n");
}

void plot_gnuplot(int mode) {
    char *cmds_gnuplot[NB_CMDS] = {};
    cmds_gnuplot[0] = "set terminal jpeg 1920,1080";
    cmds_gnuplot[2] = "set title Psi en fonction de x";
    cmds_gnuplot[3] = "set xlabel 'x'";
    cmds_gnuplot[4] = "set ylabel 'Psi'";

    switch (mode) {
    //Mode puit de potentiel infini
    case 0:
        cmds_gnuplot[1] = "set output 'plot/puit_inifini.jpeg'";
        cmds_gnuplot[5] = "plot 'data/puit_infini.dat' w l lw 1";                                 
        break;
    //
    default: break;
    }
    //Ouverture du flux vers gnuplot
    FILE *gnuplot_pipe = popen("gnuplot -persistent", "w");
    if (!gnuplot_pipe) {
        fprintf(stderr, "Error, failed to open pipe to gnuplot\n");
        exit(EXIT_FAILURE);
    }
    //Execution des commandes gnuplot
    for (int i = 0; i < NB_CMDS; i++) {

        /* PROBLEME SEGV
        fprintf(gnuplot_pipe, "%s\n", cmds_gnuplot[i]);
        */
       
    }
    pclose(gnuplot_pipe);
}

void euler_method(double x, double y[], double f[], void *params_ptr) {
    double *params = (double *)params_ptr;
    double xmax = params[0];
    double h = params[1];
    double k = params[2];
    double E = params[3];
    double V = params[4];

    while (x < xmax) {
        //ddpsi = (-k)*(E-V)*psi
        f[1] = (-k)*(E-V)*y[0];
        //dpsi = dpsi + h*ddpsi
        f[0] = f[0] + h*f[1];
        y[1] = f[0];
        //psi = psi + h*dpsi
        y[0] = y[0] + h*y[1];
        //dN = psi*psi
        f[2] = y[0]*y[0];
        //N = N + h*dN
        y[2] = y[2] + h*f[2];
        x += h;
    }
}

void solve_euler(_Psi *p, double psi, double dpsi, double m, double E, double V, double xmin, double xmax) {
    double h = xmax/NB_PTS;         //Step d'incrementation pour la resolution
    double k = 2*m/(HBARC*HBARC);   //Paramètre de l'eq diff
    double params[5] = {xmax, h, k, E, V};                 
    
    double N = 0.0;
    double y[3] = {psi, dpsi, N};

    double ddpsi = (-k)*(E-V)*psi;
    double dN = psi*psi;
    double f[3] = {dpsi, ddpsi, dN};

    euler_method(xmin, y, f, &params);
    //Si la normalisation n'est pas validee, on change dpsi pour que N se rapproche de 1
    if (fabs(y[2])-1 < 1e-4) {
        //dpsi = 1/sqrt(N)
        f[0] = (1/sqrtf(y[2]));
        y[1] = f[0];
        euler_method(xmin, y, f, &params);
    }
    // printf("psi(L) = %.12lf\ndpsi = %.12lf\nN(L) = %.12lf\nE = %lf eV\n\n**************************\n\n", y[0], y[1], y[2], E);
    p->psi = y[0];
    p->dpsi = y[1];
}

void puit_potentiel_infini(double m, double L, double E) {
    double dE = 1e-3;               //Step d'incrementation de l'energie
    double V = 0.0;                 //Potentiel à l'interieur du puit
    double xmin = 0.0, xmax = L;
    double psi = 0.0, dpsi = 1.0;
    
    double dpsi_n[5] = {0.0};       //Tableau des valeurs d'entree solutions
    double E_n[5] = {0.0};          //Tableau des niveaux d'energie quantifiables
    int count = 0; 
    
    //Recherche des 5 premieres solutions (E, dpsi) du systeme
    double psi_xmax = 1.0;
    _Psi *p = (_Psi *)calloc(1, sizeof(_Psi));
    while (count < 5) {
        while (fabs(psi_xmax) > 2e-3) {
            E += dE;
            solve_euler(p, psi, dpsi, m, E, V, xmin, xmax);
            psi_xmax = p->psi;
        }
        psi_xmax = 1.0;
        dpsi_n[count] = p->dpsi;
        E_n[count] = E;
        count++;
    }
    display_array(E_n, 5, "All 5 E");
    display_array(dpsi_n, 5, "All 5 dpsi");

    //Creation des rendus visuels avec gnuplot
    FILE *data = fopen("data/puit_infini.dat", "w");
    if (!data) {
        fprintf(stderr, "Error, failed to open %s\n", "data/puit_infini.dat");
        exit(EXIT_FAILURE);
    }
    fprintf(data, "x Psi\n");
    double h = xmax/NB_PTS;
    double k = 2*m/(HBARC*HBARC);
    // Stockage de toutes les donnees pour chaque solution (E, dpsi)
    for (int i = 0; i < 5; i++) {
        double ddpsi = (-k)*(E_n[i]-V)*psi;
        double xdata[NB_PTS] = {0.0};
        double psidata[NB_PTS] = {0.0};
        double y[2] = {psi, dpsi_n[i]};
        double f[2] = {dpsi_n[i], ddpsi};
        double params[5] = {xmax, h, k, E_n[i], V};
        store_data(xdata, psidata, xmin, y, f, &params);
        //Ecriture des toutes les donnees dans un fichier temporaire
        write_data(data, xdata, psidata);   
    }
    fclose(data);
    //Plot le fichier de donnees avec gnuplot
    plot_gnuplot(0);

    free(p);
}