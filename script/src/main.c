// ::setlocal makeprg=cd\ script\ &&\ make\ main\ &&\ ./main.x
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/TESI_fun.h"


void chack_parameters(double l, double *E, double *r0, int sign){

    printf("Running with\nl\t%.3f\n", l);
    printf("E\t%.3e\n", *E);
    printf("sign\t%d\n", sign);

    double V_data[4];
    TESI_Veff_max_min(l, V_data);

    printf("\nV_max   %.3e\n", V_data[2]);
    printf("V_min  %.3e\n", V_data[3]);
    printf("E       %.3e => ", *E);

    if (*E >= V_data[2] && sign == 1)
        printf("escape (the positive solution was chosen)\n");
    else if (*E >= V_data[2] && sign == -1)
        printf("infall\n");
    else if (*E < V_data[2] && *E > 0)
        printf("unbound orbit\n");
    else if (*E < 0 && *E > V_data[3]){
        printf("bound orbit, changing r0 to r_min\n");
        *r0 = V_data[1];
    }
    else if (*E < V_data[3])
        printf("infall and already really close!\n");
    else {
        printf("not possible\n");
        exit(1);
    }
    printf("r0\t%.3f\n", *r0);
}


void print_help(char *argv[]){
    printf("Usage: %s l E -opt1 2.5 -opt2 3...\n", argv[0]);
    printf("l:\tangular momentum per unit rest mass per Sh. radius\n");
    printf("E:\tkinetic energy per unit rest mass\n");
    printf("\nOptional arguments:\n");
    printf("-s :\tsign, direction of initial radial velocity (-1 inward, 1 outward) (default is -1)\n");
    printf("-r :\tr0, starting radius in unit of Sh. radius (10 default)\n");
    printf("-h :\th, proper time increment (1e-3 default)\n");
    printf("-t :\ttau_max, maximum proper time (100 default)\n");
    printf("\nl < sqrt(3) no stable points\n");
    printf("l = sqrt(3) one stationary point (r_ISCO = 3)\n");
    printf("l > sqrt(3) two stationary points\n");
    printf("Provide l for more information about Veff\n");
    exit(1);
}


int main(int argc, char *argv[]){

    if (argc < 2 || argc > 7)
        print_help(argv);


    double l;
    double V_data[4];
    if (argc == 2){
        l = atof(argv[1]);
        printf("l =\t\t%.3f\n", l);
        TESI_Veff_max_min(l, V_data);
        printf("r_max\t\t%.3f\n", V_data[0]);
        printf("r_min\t\t%.3f\n", V_data[1]);
        printf("V(r_max)\t%.3e\n", V_data[2]);
        printf("V(r_min)\t\b%.3e\n", V_data[3]);
        printf("Provide E to start the simulation\n");
        exit(1);
    }

    l = atof(argv[1]);
    double E = atof(argv[2]);

    // Simulation parameters
    double h = 1e-3;
    int sign = -1;
    double r0 = 10.;
    double tau_max = 100;

    double r = r0;
    double phi = 0.;
    double t = 0.;
    double tau = 0.;

    // Optional arguments check
    for (int i = 3; i < argc; i++){
        if (argv[i][0] != '-'){
            printf("Invalid argument %s\n", argv[i]);
            exit(1);
        }
        switch (argv[i][1]){
            case 's':
                sign = atoi(argv[i + 1]);
                break;
            case 'r':
                r0 = atof(argv[i + 1]);
                break;
            case 'h':
                h = atof(argv[i + 1]);
                break;
            case 't':
                tau_max = atof(argv[i + 1]);
                break;
            default:
                printf("Invalid argument %s\n", argv[i]);
                print_help(argv);
        }
        i++;
    }

    chack_parameters(l, &E, &r, sign);
    printf("tau_max\t%.3f\n", tau_max);
    printf("h\t%.3e\n\n", h);


    FILE *f = fopen("data/orbit.csv", "w");
    fprintf(f, "tau,r,phi,t\n");

    while (tau < tau_max){
        TESI_RK4(h, tau, &r, &phi, &t, E, l, &sign);
        tau += h;
        if (r < 0.1){
            printf("\nMass reached (r < 0.1). Simulation terminated.\n");
            break;
        }
        if (r > r0){
            printf("\nEscape reached (r > %.0f). Simulation terminated.\n", r0);
            break;
        }
        if (fmod(tau, 0.02) < h){
            fprintf(f, "%.10e,%.10e,%.10e,%.10e\n", tau, r, phi, t);
        }
    }

    fclose(f);

    printf("\nEnded at:\n");
    printf("r = %f\n", r);
    printf("phi = %f\n", phi / (2 * M_PI));
    printf("t = %f\n", t);
    printf("tau = %f\n", tau);
}

// Unbound orbit:   ./main.x 10 3.5
// Precession:      ./main.x 3 -0.0032 -t 15000 -r 1000 
// Infall:          ./main.x 3 0.2
