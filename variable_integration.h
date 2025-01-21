#pragma once

#include "integration.h"

//Function to calculate the total mass up to a given point of a subgrid
double mass_ir(double ***all_subgrids, int grid_n, int ir_stop){

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

    //Integrate mass from subgrids inside event horizon
    for(int ig=N_LEVELS;ig>grid_n;ig--){
        //Load each subgrid (in order)
        grid_l = all_subgrids[ig];
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
    //Load subgrid with event horizon
    grid_l = all_subgrids[grid_n];
    r = grid_l[0];
    Phi = grid_l[2];
    Pi = grid_l[3];
    a = grid_l[4];
    dR = *(grid_l[6]);
    nR = (int)(*(grid_l[8]));

    integrand_previous = 2.0*PI*(Phi[0]*Phi[0] +Pi[0]*Pi[0])*r[0]*r[0]/(a[0]*a[0]);
    for(int ir=1;ir<ir_stop+1;ir++){
        integrand = 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir]);
        mass += 0.5*dR*(integrand+integrand_previous);
        integrand_previous = integrand;
    }

    return mass;
}

//Function to integrate the metric only in a specific subgrid
void partial_metric(double ***all_subgrids, int n_levels, int grid_n){
    //Calculate the full metric (without normalization)
    //Load the current grid
    double **grid_l = all_subgrids[grid_n];
    double *r = grid_l[0];
    double *Phi = grid_l[2];
    double *Pi  = grid_l[3];
    double *a     = grid_l[4];
    double *alpha = grid_l[5];
    double deltaR =   *(grid_l[6]);
    int nR = (int)(*(grid_l[8]));
    
    //Calculate the sum of the metric for the previous subgrids
    double     a0 = 1.0;
    double alpha0 = 1.0;
    for (int ig=n_levels; ig>grid_n; ig--){
        //Load the partial sums of the subgrid
        double partial_a     = all_subgrids[ig][13][0];
        double partial_alpha = all_subgrids[ig][14][0];
            a0 += partial_a;
        alpha0 += partial_alpha;
    }

    //Calculate auxiliar variable Beta
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
    metric_iteration(a0, alpha0, Beta, Beta_p05, a, alpha, r, nR, deltaR, 1.0);
    
    //Calculate the sum of the metric for the following subgrids
    double     a_final =     a[nR-EXTRA_N];
    double alpha_final = alpha[nR-EXTRA_N];
    for (int ig=grid_n-1; ig>0; ig--){
        //Load the partial sums of the subgrid
        double partial_a     = all_subgrids[ig][13][0];
        double partial_alpha = all_subgrids[ig][14][0];
            a_final += partial_a;
        alpha_final += partial_alpha;
    }

    //Normalize the metric
    double norm = 1.0/(a_final*alpha_final);
    for(int ir=0;ir<nR;ir++){
        alpha[ir] = alpha[ir]*norm;
    }
}

//Function to integrate the whole metric
void full_metric(int save_levels,int n_levels, double ***all_subgrids){
    //printf("call to metric\n");
    double a0 = 1.0;
    double alpha0 = 1.0;
    double *a;
    double *alpha;
    double deltaR;
    int nR;
    //Calculate the full metric (without normalization)
    for (int i=0; i<n_levels; i++){

        //Load the current grid
        double **grid_l = all_subgrids[n_levels-i];
        double *r = grid_l[0];
        double *Phi = grid_l[2];
        double *Pi  = grid_l[3];
        a     = grid_l[4];
        alpha = grid_l[5];
        deltaR =   *(grid_l[6]);
        nR = (int)(*(grid_l[8]));
        double *partial_a_p = grid_l[13];
        double *partial_alpha_p = grid_l[14];
        //printf("grid %d; deltaR %lf; nR %d\n",n_levels-i,deltaR,nR);

        //Calculate auxiliar variable Beta
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
        metric_iteration(a0, alpha0, Beta, Beta_p05, a, alpha, r, nR, deltaR, 1.0);
        if(i == n_levels-1){
            a0     = a[nR-1];
            alpha0 = alpha[nR-1];
        } else {
            a0     = a[nR-EXTRA_N];
            alpha0 = alpha[nR-EXTRA_N];
        }
        partial_a_p[0] = a0 - a[0];
        partial_alpha_p[0] = alpha0 - alpha[0];
        //printf("a0: %e; alpha0 %e\n\n",a0,alpha0);
    }
    //Normalize the metric
    double norm = 1.0/(a0*alpha0);
    norm_global = norm;
    for (int i=1;i<n_levels+1;i++){
        alpha = all_subgrids[i][5];
        nR = (int)(*(all_subgrids[i][8]));
        for(int ir=0;ir<nR;ir++){
            alpha[ir] = alpha[ir]*norm;
        }
    }
}