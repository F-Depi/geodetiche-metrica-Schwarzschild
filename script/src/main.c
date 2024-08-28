// ::setlocal makeprg=cd\ script\ &&\ make\ main\ &&\ ./main.x
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/TESI_fun.h"


void chack_parameters(double l, double E, int sign){

    printf("Running with l = %.3f and E = %.3f\n", l, E);

    double V_data[4];
    TESI_Veff_max_min(l, V_data);

    printf("\nV_max   %.3f\n", V_data[2]);
    printf("V_min  %.3f\n", V_data[3]);
    printf("E       %.3f => ", E);

    if (E >= V_data[2] && sign == 1)
        printf("escape (the positive solution was chosen)\n");
    else if (E >= V_data[2] && sign == -1)
        printf("infall\n");
    else if (E < V_data[2] && E > 0)
        printf("unbound orbit\n");
    else if (E < 0 && E > V_data[3])
        printf("bound orbit\n");
    else if (E < V_data[3])
        printf("infall and already really close!\n");
    else {
        printf("not possible\n");
        exit(1);
    }
}


int main(int argc, char *argv[]){

    if (argc > 4 || argc < 3){
        printf("Usage: %s l E\n", argv[0]);
        printf("r0 is optional, 10 (or less is necessary is default)\n");
        exit(1);
    }
    double l = atof(argv[1]);
    double E = atof(argv[2]);
    double r0 = 10.;
    double r = r0;
    double phi = 0.;
    double t = 0.;
    double tau = 0.;
    double h = 1e-3;
    int sign = -1;

    chack_parameters(l, E, sign);

    FILE *f = fopen("data/orbit.csv", "w");
    fprintf(f, "tau,r,phi,t\n");

    while (r > 0.1 && r <= 1.1 * r0 && (tau / h) < 1e5){
        TESI_RK4(h, tau, &r, &phi, &t, E, l, &sign);
        tau += h;
        fprintf(f, "%.10e,%.10e,%.10e,%.10e\n", tau, r, phi, t);
    }

    fclose(f);

    printf("\nEnded at:\n");
    printf("r = %f\n", r);
    printf("phi = %f\n", phi / (2 * M_PI));
    printf("t = %f\n", t);
    printf("tau = %f\n", tau);
}

