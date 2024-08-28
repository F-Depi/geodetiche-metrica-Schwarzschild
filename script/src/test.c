// ::setlocal makeprg=cd\ script\ &&\ make\ test\ &&\ ./test.x
#include <stdio.h>
#include <math.h>
#include "../include/TESI_fun.h"


void t_hard_to_compute_well(){
    double l = 0;
    double E = 0;

    double r = 1.1;
    double phi = 0.;
    double t = 0.;
    double tau = 0.;
    double h = 1e-8;
    int sign = -1;

    FILE *f = fopen("data/test_infall.csv", "w");
    fprintf(f, "tau,r,phi,t\n");

    while (r > 1.){
        TESI_RK4(h, tau, &r, &phi, &t, E, l, &sign);
        tau += h;
        if (fmod(tau, 0.001) < 1e-9) 
            fprintf(f, "%.3f,%.3f,%.3f,%.3f\n", tau, r, phi, t);
    }

    fclose(f);
    printf("r = %f\n", r);
    printf("phi = %f\n", phi);
    printf("t = %f\n", t);
    printf("tau = %f\n", tau);
}


int main(){

    double l = sqrt(3);
    double r_ISCO = 3;

    printf("%f\n", TESI_Veff(r_ISCO, l));

    FILE *f = fopen("data/Veff.csv", "w");
    fprintf(f, "r,Veff\n");

    for (double r = 1; r < 10; r += 0.01){
        fprintf(f, "%.6e,%.6e\n", r, TESI_Veff(r, l));
    }

    fclose(f);


    return 0;
}
