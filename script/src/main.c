// ::setlocal makeprg=cd\ script\ &&\ make\ main\ &&\ ./main.x
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
    printf("V(r_min)\t%.6e\n", V_data[3]);
    printf("Provide E to start the simulation\n");
    exit(1);
}


// Check E from argv[2], if it is "max" or "min" calculate the local maximum or
// minimum of Veff and assign it to E
double parse_E(char *argv[]){
    double V_data[4];
    TESI_Veff_max_min(atof(argv[1]), V_data);
    
    if (strcmp(argv[2], "max") == 0){
        return V_data[2];
    }
    else if (strcmp(argv[2], "min") == 0){
        return V_data[3];
    }
    return atof(argv[2]);
}


// Write the filename for the output data, using l and E from argv
// If E is "max" or "min" the filename will be l%.3f_Emax.csv or l%.3f_Emin.csv
void write_filename(char *filename, char *argv[]){
    if (strcmp(argv[2], "max") == 0){
        sprintf(filename, "data/l%.3f_Emax.csv", atof(argv[1]));
        return;
    }
    else if (strcmp(argv[2], "min") == 0){
        sprintf(filename, "data/l%.3f_Emin.csv", atof(argv[1]));
        return;
    }
    sprintf(filename, "data/l%.3f_E%.5f.csv", atof(argv[1]), atof(argv[2]));
}


int main(int argc, char *argv[]){

    if (argc < 2){
        printf("Missing arguments\n");
        print_help(argv);
    }
    if (argc > 15){
        printf("Too many arguments\n");
        print_help(argv);
    }

    if (argc == 2)
        print_help2(argv);


    /***** Orbit parameters *****/
    double l = atof(argv[1]);       // Angular momentum per unit rest mass per Sh. radius
    double E = parse_E(argv);       // Kinetic energy per unit rest mass
    double V_data[4];               // Veff(r) r max and min + V values
    double r0 = 10.;                // Starting radius in unit of Sh. radius
    double r_lim = 10.;             // Biggest r for the simulation
    int sign = -1;                  // Initial radial velocity direction
    int Nturns = 0;                 // Number of turns

                                    
    /***** Other parameters that can be changed with optional arguments *****/
    double h = 1e-3;                // Proper time increment
    double tau_max = 100;           // Maximum proper time
    int time2print;                 // Data saved every time2print steps
    char filename[50];              // Output file
    write_filename(filename, argv);
    char fps[20];                   // Frames per second
    sprintf(fps, "50");

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
            case 'f':
                sprintf(filename, "data/%s.csv", argv[i + 1]);
                break;
            case 'B':
                sprintf(fps, "%s", argv[i + 1]);
                break;
            default:
                printf("Invalid argument %s\n", argv[i]);
                print_help(argv);
        }
        i++;
    }

    if (strcmp(fps, "all") == 0) 
        time2print = 1;
    else 
        time2print = ceil(1 / h / atof(fps));

    check_parameters(l, E, &r0, &r_lim, &sign);
    printf("r_lim\t%.3f\n", r_lim);
    printf("h\t%.3e\n", h);
    printf("tau_max\t%.3f\n", tau_max);
    printf("Saving data every %d steps (fps = %s)\n", time2print, fps);
    printf("Output file: %s\n\n", filename);


    /***** Coordinates *****/
    double tau = 0.;                // Proper time
    double r = r0;                  // Radius
    double phi = 0.;                // Azimuthal angle
    double t = 0.;                  // Schwarzschild time

    
    FILE *f = fopen(filename, "w");
    fprintf(f, "tau,r,phi,t\n");
    fprintf(f, "%.10e,%.10e,%.10e,%.10e\n", tau, r, phi, t);

    int kk = 0;
    while (tau < tau_max){

        TESI_RK4(h, tau, &r, &phi, &t, E, l, &sign, &Nturns);
        tau += h;
        kk++;

        if (kk % time2print == 0){
            fprintf(f, "%.15e,%.10e,%.10e,%.15e\n", tau, r, phi, t);
            printf("\rtau = %.3e | r = %.3f | Turns = %d", tau, r, Nturns);
            fflush(f);
        }

        if (r < 0.5){
            printf("\nMass reached (r < 0.1). Simulation terminated.\n");
            break;
        }
        if (r > 1.1 * r_lim){
            printf("\nEscape reached (r > %.0f). Simulation terminated.\n", r_lim);
            break;
        }
    }

    fclose(f);

    printf("\nEnded at:\n");
    printf("tau = %f\n", tau);
    printf("r = %f\n", r);
    printf("phi = %f\n", phi / (2 * M_PI));
    printf("t = %f\n", t);

    return 0;
}
