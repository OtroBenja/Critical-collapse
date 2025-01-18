#pragma once

#include "integration.h"

//Function to calculate the total mass up to a given point of a subgrid
double mass_ir(double ***all_subgrids, int grid_n, int ir);

//Function to integrate the whole metric
void full_metric(int level,int n_levels, double ***all_subgrids);



//Function to calculate the total mass up to a given point of a subgrid
double mass_ir(double ***all_subgrids, int grid_n, int ir){

    double **grid_l;
    double *r;
    double *Phi;
    double *Pi;
    double *a;
    double dR;
    int nR;

    double integrand;
    double integrand_previous;
    double temp;
    double mass = 0;

    for(int n=0;n<grid_n;n++){
        //Load each subgrid (in order)
        grid_l = all_subgrids[N_LEVELS-n];
        r = grid_l[0];
        Phi = grid_l[2];
        Pi = grid_l[3];
        a = grid_l[4];
        dR = *(grid_l[6]);
        nR = (int)(*(grid_l[8]));

        integrand_previous = 2.0*PI*(Phi[0]*Phi[0] +Pi[0]*Pi[0])*r[0]*r[0]/(a[0]*a[0]);
        for(int ir=1;ir<nR-EXTRA_N+1;ir++){
            integrand = 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir]);
            mass += 0.5*dR*(integrand+integrand_previous);
            integrand_previous = integrand;
        }
    }
    //Load each subgrid (in order)
    grid_l = all_subgrids[N_LEVELS-grid_n];
    r = grid_l[0];
    Phi = grid_l[2];
    Pi = grid_l[3];
    a = grid_l[4];
    dR = *(grid_l[6]);
    nR = (int)(*(grid_l[8]));

    integrand_previous = 2.0*PI*(Phi[0]*Phi[0] +Pi[0]*Pi[0])*r[0]*r[0]/(a[0]*a[0]);
    for(int ir=1;ir<nR+1;ir++){
        integrand = 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir]);
        mass += 0.5*dR*(integrand+integrand_previous);
        integrand_previous = integrand;
    }

    return mass;
}

//Function to integrate the whole metric
void full_metric(int save_levels,int n_levels, double ***all_subgrids){
    //printf("call to metric\n");
    double a0 = 1.0;
    double alpha0 = 1.0;
    double *r;
    double *Phi;
    double *Pi;
    double *a;
    double *alpha;
    double deltaR;
    int nR;
    //Calculate the full metric (without normalization)
    for (int i=0; i<n_levels; i++){

        //Load the current grid
        double **grid_l = all_subgrids[n_levels-i];
        r = grid_l[0];
        Phi = grid_l[2];
        Pi  = grid_l[3];
        a     = grid_l[4];
        alpha = grid_l[5];
        deltaR =   *(grid_l[6]);
        nR = (int)(*(grid_l[8]));
        //printf("grid %d; deltaR %lf; nR %d\n",n_levels-i,deltaR,nR);

        //Calculate auxiliar variable Beta
        //double     *Beta = malloc(sizeof(double)*nR);
        //double *Beta_p05 = malloc(sizeof(double)*nR);
        double     Beta[nR];
        double Beta_p05[nR];
        for(int ir=0;ir<nR-1;ir++){
            Beta[ir] = 2.0*PI*r[ir]*((Pi[ir])*(Pi[ir]) + (Phi[ir])*(Phi[ir]));
            Beta_p05[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*(
                ( Pi[ir] + Pi[ir+1])*( Pi[ir] + Pi[ir+1]) +
                (Phi[ir] +Phi[ir+1])*(Phi[ir] +Phi[ir+1]));
        }
        Beta[nR-1] = 2.0*PI*r[nR-1]*((Pi[nR-1])*(Pi[nR-1]) + (Phi[nR-1])*(Phi[nR-1]));
        //Iterate metric and save final value to continue on the next grid
        //for(int j=0;j<10;j++){
        //    printf("r%d: %e\n",j,r[j]);
        //    printf("Phi%d: %lf; Pi%d %lf\n" ,j,Phi[j],j,Pi[j]);
        //    printf("a%d: %lf; alpha%d %lf\n",j,  a[j],j,alpha[j]);
        //}
        //if(i<save_levels){
        metric_iteration(a0, alpha0, Beta, Beta_p05, a, alpha, r, nR, deltaR, 1.0);
        a0     = a[nR-EXTRA_N];
        alpha0 = alpha[nR-EXTRA_N];
        //} else {metric_iteration_NW(a0, alpha0, Beta, Beta_p05, &a0, &alpha0, r, nR+1-EXTRA_N, deltaR, 1.0);}
        //printf("a0: %e; alpha0 %e\n\n",a0,alpha0);

        //free(Beta);
        //free(Beta_p05);
    }
    //Normalize the metric
    double norm = 1.0/(a0*alpha0);
    for (int i=1;i<n_levels+1;i++){
        alpha = all_subgrids[i][5];
        nR = (int)(*(all_subgrids[i][8]));
        for(int ir=0;ir<nR;ir++){
            alpha[ir] = alpha[ir]*norm;
        }
    }
}