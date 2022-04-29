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

/********** Solving Time-Independent Schrodinger equation **********/

void store_data(double xdata[], double psidata[], double x, double y[], double f[], void *params_ptr) {
    double *params = (double *)params_ptr;
    double xmax = params[0];
    double h = params[1];
    double k = params[2];
    double E = params[3];
    double V = params[4];

    psidata[0] = y[0];
    xdata[0] = x;
    int i = 1;
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
        // fprintf(data, "%lf %lf\n", xdata[i], fabs(psidata[i]));
        fprintf(data, "%lf %lf\n", xdata[i], psidata[i]);
    }
    fprintf(data, "\n");
}

void plot_gnuplot(int mode) {
    char *cmds_gnuplot[NB_CMDS] = {};
    cmds_gnuplot[0] = "set terminal jpeg size 1920,1080";
    cmds_gnuplot[2] = "set title \"Psi en fonction de x\"";
    cmds_gnuplot[3] = "set xlabel \"x\"";
    cmds_gnuplot[4] = "set ylabel \"Psi\"";

    switch (mode) {
    case 0:     //Infinite potential well
        cmds_gnuplot[1] = "set output \"plot/puit_inifini.jpeg\"";
        cmds_gnuplot[5] = "plot \"data/puit_infini.dat\" w l title \"Psi(x)\"";                                 
        break;
    case 1:     //
        break;
    case 2:     //
        break;
    default: break;
    }
    FILE *setup_gnuplot = fopen("data/setup.dat", "w");
    if (!setup_gnuplot) {
        fprintf(stderr, "Error, failed to open setup_gnuplot\n");
        exit(EXIT_FAILURE);
    }
    //Writing gnuplot commands un setup file
    for (int i = 0; i < NB_CMDS; i++) {
        fprintf(setup_gnuplot, "%s\n", cmds_gnuplot[i]);
    }
    fclose(setup_gnuplot);
    system("gnuplot -p < data/setup.dat");
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
    double h = xmax/NB_PTS;         //Increment step for solving the system
    double k = 2*m/(HBARC*HBARC);   //Parameter of the differential equation
    double params[5] = {xmax, h, k, E, V};                 
    
    double N = 0.0;
    double y[3] = {psi, dpsi, N};

    double ddpsi = (-k)*(E-V)*psi;
    double dN = psi*psi;
    double f[3] = {dpsi, ddpsi, dN};

    euler_method(xmin, y, f, &params);
    //If normalisation isn't checked, change dspi to make N close to 1
    if (fabs(y[2])-1 < 1e-4) {
        //dpsi = 1/sqrt(N)
        f[0] = fabs(1/sqrtf(y[2]));
        y[1] = f[0];
        euler_method(xmin, y, f, &params);
    }
    p->psi = y[0];
    p->dpsi = y[1];
}

void infinite_potential_well(double m, double L, double E) {
    double dE = 1e-4;               //Increment step for energy
    double V = 0.0;                 //Potential inside the well
    double xmin = 0.0, xmax = L;
    double psi = 0.0, dpsi = 1.0;
    
    double dpsi_n[5] = {0.0};       //Array of solutions dpsi
    double E_n[5] = {0.0};          //Array of quantifiable energy levels
    int count = 0; 
    
    //Search 5 first solutions (E, dpsi) of the system
    double psi_xmax = 1.0;
    _Psi *p = (_Psi *)calloc(1, sizeof(_Psi));
    while (count < 5) {
        while (fabs(psi_xmax) > 2e-4) {
            E += dE;
            solve_euler(p, psi, dpsi, m, E, V, xmin, xmax);
            psi_xmax = p->psi;
        }
        psi_xmax = 1.0;
        dpsi_n[count] = fabs(p->dpsi);
        E_n[count] = E;
        count++;
    }
    display_array(E_n, 5, "All 5 E");
    display_array(dpsi_n, 5, "All 5 dpsi");
    //Generating graphs using gnuplot
    FILE *data = fopen("data/puit_infini.dat", "w");
    if (!data) {
        fprintf(stderr, "Error, failed to open %s\n", "data/puit_infini.dat");
        exit(EXIT_FAILURE);
    }
    double h = xmax/NB_PTS;
    double k = 2*m/(HBARC*HBARC);
    //Store all data for each solution (E, dpsi)
    for (int i = 0; i < 5; i++) {
        double ddpsi = (-k)*(E_n[i]-V)*psi;
        double xdata[NB_PTS] = {0.0};
        double psidata[NB_PTS] = {0.0};
        double y[2] = {psi, dpsi_n[i]};
        double f[2] = {dpsi_n[i], ddpsi};
        double params[5] = {xmax, h, k, E_n[i], V};
        store_data(xdata, psidata, xmin, y, f, &params);
        //Write all data in temporary file
        write_data(data, xdata, psidata);   
    }
    fclose(data);
    plot_gnuplot(0);
    printf("Successfully plot to file [plot/puit_infini.jpeg]\n");

    free(p);
}