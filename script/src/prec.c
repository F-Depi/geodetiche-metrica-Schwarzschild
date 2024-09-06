// ::setlocal makeprg=cd\ script\ &&\ make\ prec\ &&\ ./prec.x
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../include/TESI_fun.h"


void print_help(char *argv[]){
    printf("Usage: %s l E -opt1 2.5 -opt2 3...\n", argv[0]);
    printf("l:\tangular momentum per unit rest mass per Sh. radius\n");
    printf("E:\tkinetic energy per unit rest mass\n");
    printf("\nOptional arguments:\n");
    printf("-s :\tsign, direction of initial radial velocity (-1 inward, 1 outward) (default is -1)\n");
    printf("-r :\tr0, starting radius in unit of Sh. radius (10 default)\n");
    printf("-h :\th, proper time increment (1e-3 default)\n");
    printf("-t :\ttau_max, maximum proper time (100 default)\n");
    printf("-f :\tfilename, output file (default is data/l\%.3f_E\%.5f.csv default)\n");
    printf("-B :\tframes per second (50 default, all to save all data)\n");
    printf("\nl < sqrt(3) no stable points\n");
    printf("l = sqrt(3) one stationary point (r_ISCO = 3)\n");
    printf("l > sqrt(3) two stationary points\n");
    printf("Provide l for more information about Veff\n");
    exit(1);
}


void print_help2(char *argv[]){
    double l = atof(argv[1]);
    double V_data[4];
    printf("l =\t\t%.6f\n", l);
    TESI_Veff_max_min(l, V_data);
    printf("r_max\t\t%.6f\n", V_data[0]);
    printf("r_min\t\t%.6f\n", V_data[1]);
    printf("V(r_max)\t%.6e\n", V_data[2]);
    printf("V(r_min)\t\b%.6e\n", V_data[3]);
    printf("Provide E to start the simulation\n");
    exit(1);
}


int main(int argc, char *argv[]){

    if (argc < 2){
        printf("Missing arguments\n");
        print_help(argv);
    }
    if (argc > 11){
        printf("Too many arguments\n");
        print_help(argv);
    }

    if (argc == 2)
        print_help2(argv);


    /***** Orbit parameters *****/
    double l = atof(argv[1]);       // Angular momentum per unit rest mass per Sh. radius
    double E = atof(argv[2]);       // Kinetic energy per unit rest mass
    double V_data[4];               // Veff(r) r max and min + V values
    double r0 = 10.;                // Starting radius in unit of Sh. radius
    double r_lim = 10.;             // Biggest r for the simulation
    int sign = -1;                  // Initial radial velocity direction
    int Nturns = 0;                 // Number of turns

                                    
    /***** Other parameters that can be changed with optional arguments *****/
    double h = 1e-3;                // Proper time increment
    double tau_max = 100;           // Maximum proper time

    /***** Optional arguments check *****/
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
                r_lim = atof(argv[i + 1]);
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
    
    check_parameters(l, E, &r0, &r_lim, &sign);
    printf("h\t%.3e\n", h);
    printf("tau_max\t%.3f\n\n", tau_max);


    /***** Coordinates *****/
    double tau = 0.;                // Proper time
    double r = r0;                  // Radius
    double phi = 0.;                // Azimuthal angle
    double t = 0.;                  // Schwarzschild time

    
    /***** Precession *****/
    char filename_prec[50];
    sprintf(filename_prec, "data/precession/l%.3f_E%.5f_h%.0e.csv", l, E, h);
    FILE *f_prec = fopen(filename_prec, "w");
    fprintf(f_prec, "r,phi\n");
    int Nturns_old;
    double r_old, r_old_old;
    double phi_old;

    int kk = 0;

    double v = TESI_fun_r(r, E, l, &sign, &Nturns);
    double E0 = v*v/2. + TESI_Veff(r, l);

    while (tau < tau_max){

        Nturns_old = Nturns;
        r_old_old = r_old;
        r_old = r;
        phi_old = phi;

        // TESI_RK4(h, tau, &r, &phi, &t, E, l, &sign, &Nturns);
        TESI_RK4_corrected(h, tau, &v, &r, &phi, &t, E, l, &sign, &Nturns);
        // TESI_RK4_corrected2(h, tau, &r, &phi, &t, E, l, &sign, &Nturns);
        // TESI_RKN4(h, tau, &v, &r, &phi, &t, E, l, &sign, &Nturns);
        tau += h;
        kk++;

        if (r < 1){
            printf("\nMass reached (r < 0.1). Simulation terminated.\n");
            break;
        }

        if (r_old_old <= r_old && r_old >= r){
            printf("\nOuter turning point at r = %.3f\n", r_old);
            fprintf(f_prec, "%.15e,%.15e\n", r_old, phi_old);
        }

        if (r_old_old >= r_old && r_old <= r){
            printf("\nInner turning point at r = %.3f\n", r_old);
            fprintf(f_prec, "%.15e,%.15e\n", r_old, phi_old);
        }

        if (kk % 100 == 0)
            printf("\rtau = %.3e | r = %.3f | Turns = %d", tau, r, Nturns);

        if (r > 1.1*r_lim){
            printf("\nEscape reached (r > %.0f). Simulation terminated.\n", r_lim);
            break;
        }
    }

    if (tau >= tau_max)
        printf("\nMaximum proper time reached\n");

    double E_fin = v*v/2. + TESI_Veff(r, l);
    printf("E = %.7e\n", E_fin);
    printf("E - E0 = %.7e\n", (E_fin - E0) / E0);
    fclose(f_prec);

    return 0;
}
