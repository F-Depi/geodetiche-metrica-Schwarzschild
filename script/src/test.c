// ::setlocal makeprg=cd\ script\ &&\ make\ test\ &&\ ./test.x
#include <stdio.h>
#include <math.h>
#include "../include/TESI_fun.h"


void test_Veff(double l){
    FILE *f = fopen("data/Veff.csv", "w");
    fprintf(f, "r,V\n");
    for (double r = 1; r < 100; r += 0.01){
        double V = TESI_Veff(r, l);
        fprintf(f, "%.10e,%.10e\n", r, V);
    }
    fclose(f);
}

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
        TESI_RK4(h, tau, &r, &phi, &t, E, l, &sign, NULL);
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


void test_turning_points(){
    double l = 4;
    double E = -0.003;
    double r12[2];
    TESI_turning_points(l, E, r12);
    printf("r1 = %f\n", r12[0]);
    printf("r2 = %f\n", r12[1]);
    printf("dr/dtau(r1)= %f\n", E - TESI_Veff(r12[0], l));
    printf("dr/dtau(r2) = %f\n", E - TESI_Veff(r12[1], l));
}


int main(){

    test_Veff(3);

    // t_hard_to_compute_well();

    // test_turning_points();

    return 0;
}
