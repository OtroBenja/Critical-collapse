#pragma once

#include "constants.h"
#include "integration.h"
#include "variable_integration.h"
//#include "derivatives.h"
#include "subgrid_iteration.h"
#include "load_ghostzones.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


double **initialize_subgrid(double,double*,double,double,double);
void check_collapse(double ***all_subgrids, double max_a, int iteration);


float **variable_iteration(int fType,double *model_params,double deltaR,double maxR,int iterations,int save_iteration){
    int noSaveNR = MIN_R/deltaR;
    int nR = (int)(maxR/deltaR);
    int saveNR = (int)((fmin(maxR,MAX_R)-MIN_R)/deltaR);
    double deltaT = deltaR/5.;
    double *mass;
    double max_a = pow(TOLERANCE,-0.5); // Max value of 'a' according to g^rr tolerance
    bh_mass = 0;
    bh_radius = 0;
    is_bh = false;
    last_iteration = iterations;

    //Initialize main grid
    double ***all_subgrids = malloc(sizeof(double**)*10);
    double **grid0;
    grid0 = initialize_subgrid(fType,model_params,0.0,maxR+deltaR*((double)EXTRA_N),deltaR);
    all_subgrids[0] = grid0;


    //Initialize all subgrids
    int n_levels = N_LEVELS;
    for(int n=1; n<n_levels+1; n++){
        double **gridN;
        double minR_l;
        double maxR_l;
        
        if(n==n_levels){minR_l = 0;}
        else{minR_l = pow(GRID_RATIO,-n)*maxR;}

        if(n==1){maxR_l = maxR+deltaR*((double)EXTRA_N);}
        else{maxR_l = pow(GRID_RATIO,-(n-1))*(maxR+deltaR*((double)EXTRA_N));}
        double deltaR_l = pow(GRID_RATIO,-(n-1))*deltaR;
        gridN = initialize_subgrid(fType,model_params,minR_l,maxR_l,deltaR_l);
        all_subgrids[n] = gridN;
    }
    //Initialize all subgrids

    
    //int n_levels = N_LEVELS;
    //double **gridN;
    //double deltaR_l = deltaR;
    //double maxR_l = maxR+deltaR;
    //double minR_l = 46;
    //gridN = initialize_subgrid(fType,model_params,minR_l,maxR_l,deltaR_l);
    //all_subgrids[1] = gridN;
    //deltaR_l = deltaR*0.5;
    //maxR_l = minR_l + deltaR_l;
    //minR_l = 0;
    //gridN = initialize_subgrid(fType,model_params,minR_l,maxR_l,deltaR_l);
    //all_subgrids[2] = gridN;

    double temp;
    double *temp_pointer;
    int history_size;
    float *Rhistory, *Fhistory, *Xhistory, *Yhistory, *Ahistory, *Bhistory, *Mhistory;
    if(SAVE_MODE == 0){
        Rhistory = malloc(sizeof(float)*(nR/SAVE_RES));
        history_size = sizeof(float)*(nR/SAVE_RES)*(iterations/save_iteration);
    }
    else if(SAVE_MODE == 1){
        Rhistory = malloc(sizeof(float)*saveNR);
        history_size = sizeof(float)*saveNR*(iterations-FIRST_ITERATION);
    }
    Fhistory = malloc(history_size);
    Xhistory = malloc(history_size);
    Yhistory = malloc(history_size);
    Ahistory = malloc(history_size);
    Bhistory = malloc(history_size);
    Mhistory = malloc(history_size);
    float** hist = malloc(sizeof(float*)*7);
    int save_count = save_iteration;

    //Calculate the metric before starting the iteration
    double *Beta = malloc(sizeof(double)*nR);
    double *Beta1_2 = malloc(sizeof(double)*nR);
    double *r, *phi, *Phi, *Pi, *a, *alpha;
    r     = grid0[0];
    phi   = grid0[1];
    Phi   = grid0[2];
    Pi    = grid0[3];
    a     = grid0[4];
    alpha = grid0[5];
    for(int ir=0;ir<nR;ir++){
        Beta[ir] = 2.0*PI*r[ir]*((Pi[ir])*(Pi[ir]) + (Phi[ir])*(Phi[ir]));
    }
    for(int ir=0;ir<nR-1;ir++){
        Beta1_2[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*(
            ( Pi[ir] + Pi[ir+1])*( Pi[ir] + Pi[ir+1]) +
            (Phi[ir] +Phi[ir+1])*(Phi[ir] +Phi[ir+1]));
    }
    metric_iteration(1.0,1.0,Beta, Beta1_2, a, alpha, r, nR, deltaR, 0);
    free(Beta);
    free(Beta1_2);
    full_metric(N_LEVELS,N_LEVELS,all_subgrids);

    printf("iteration start\n");
    for(int i=0;i<iterations;i++){

        if(SAVE_MODE == 0){
            //Save values of Phi, Pi, phi, a and alpha
            if(save_count == save_iteration){
                //Copy values to grid L0
                int ir = 0;
                for(int i=0; i<n_levels; i++){
                    double** gridN = all_subgrids[n_levels-i];
                    double* r_l   = gridN[0];
                    double* phi_l   = gridN[1];
                    double* Phi_l   = gridN[2];
                    double* Pi_l    = gridN[3];
                    double* a_l     = gridN[4];
                    double* alpha_l = gridN[5];
                    int nR_l = (int)(*gridN[8]);
                    for(int ir_l=0; ir_l<nR_l-EXTRA_N; ir_l += pow(GRID_RATIO,n_levels-(i+1))){
                        phi[ir]   = phi_l[ir_l];
                        Phi[ir]   = Phi_l[ir_l];
                        Pi[ir]    = Pi_l[ir_l];
                        Pi[ir]    = Pi_l[ir_l];
                        a[ir]     = a_l[ir_l];
                        alpha[ir] = alpha_l[ir_l];
                        ir++;
                    }
                }
                mass = get_mass(r,Phi,Pi,a,maxR,deltaR);
                for(int ir=0;ir<(nR/SAVE_RES);ir++){
                    Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = (float)  Phi[ir*SAVE_RES];
                    Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = (float)   Pi[ir*SAVE_RES];
                    Fhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = (float)  phi[ir*SAVE_RES];
                    Ahistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = (float)    a[ir*SAVE_RES];
                    Bhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = (float)alpha[ir*SAVE_RES];
                    Mhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = (float) mass[ir*SAVE_RES];
                }
                free(mass);
                save_count=0;
            }
            save_count+=1;
        }

        if(SAVE_MODE == 1){
            //Save values of Phi, Pi, phi, a and alpha
            if(i >= FIRST_ITERATION){
                //Copy values to grid L0
                int ir = 0;
                for(int i=0; i<n_levels; i++){
                    double** gridN = all_subgrids[n_levels-i];
                    double* phi_l   = gridN[1];
                    double* Phi_l   = gridN[2];
                    double* Pi_l    = gridN[3];
                    double* a_l     = gridN[4];
                    double* alpha_l = gridN[5];
                    int nR_l = (int)(*gridN[8]);
                    for(int ir_l=0; ir_l<nR_l-EXTRA_N; ir_l += pow(GRID_RATIO,n_levels-(i+1))){
                        phi[ir]   = phi_l[ir_l];
                        Phi[ir]   = Phi_l[ir_l];
                        Pi[ir]    = Pi_l[ir_l];
                        a[ir]     = a_l[ir_l];
                        alpha[ir] = alpha_l[ir_l];
                        ir++;
                    }
                }
                int save_iter = i - FIRST_ITERATION;
                //printf("iteration %d\n",i);
                mass = get_mass(r,Phi,Pi,a,maxR,deltaR);
                for(int ir=0;ir<saveNR;ir++){
                    Xhistory[save_iter*saveNR + ir] = (float)  Phi[noSaveNR + ir];
                    Yhistory[save_iter*saveNR + ir] = (float)   Pi[noSaveNR + ir];
                    Fhistory[save_iter*saveNR + ir] = (float)  phi[noSaveNR + ir];
                    Ahistory[save_iter*saveNR + ir] = (float)    a[noSaveNR + ir];
                    Bhistory[save_iter*saveNR + ir] = (float)alpha[noSaveNR + ir];
                    Mhistory[save_iter*saveNR + ir] = (float) mass[noSaveNR + ir];
                }
                free(mass);
            }
        }

        //Check if a collapse has happened if a > 1/TOLERANCE^2
        check_collapse(all_subgrids, max_a, i);
        //for(int ir=0;ir<nR;ir++){
        //    if(a[ir]>max_a){
        //        bh_radius = r[ir];
        //        mass = get_mass(r,Phi,Pi,a,maxR,deltaR);
        //        bh_mass = mass[ir];
        //        last_iteration = i;
        //        is_bh = true;
        //        break;
        //    }
        //}
        if(is_bh) break;
        //printf("print2\n");
        //Integrate from smallest to bigger grid
        int count_level[n_levels-1];
        for(int i=0;i<n_levels;i++){
            count_level[i] = 0;
        }
        
        for(int i=0;i<n_levels-1;i++){
            avg_overlap(all_subgrids,n_levels-i);
            left_ghostzones(all_subgrids,n_levels-i-1);
            right_ghostzones(all_subgrids,n_levels-i);
        }
        for(int i=0;i<n_levels-1;){
            //printf("level %d\n",i);
            subgrid_iteration(all_subgrids,n_levels-i,1,1,count_level[i]);
            count_level[i]++;
            if(count_level[i]==GRID_RATIO){
                count_level[i]=0;
                i++;
            } else {
                //Load previous ghostzones
                for(int j=0;j<i;j++){
                    avg_overlap(all_subgrids,n_levels-j);
                    left_ghostzones(all_subgrids,n_levels-(j+1));
                    right_ghostzones(all_subgrids,n_levels-j);
                }
                i = 0;
            }
        }
        subgrid_iteration(all_subgrids,1,1,1,1);
        full_metric(N_LEVELS,N_LEVELS,all_subgrids);

        //printf("iteration %d complete\n",i);
    }
    if(SAVE_MODE == 0){
        for(int ir=0;ir<(nR/SAVE_RES);ir++)
            Rhistory[ir] = r[ir*SAVE_RES];
    }
    else if(SAVE_MODE == 1){
        for(int ir=0;ir<(saveNR);ir++)
            Rhistory[ir] = r[noSaveNR+ir];
    }
    hist[0] = Rhistory;
    hist[1] = Fhistory;
    hist[2] = Xhistory;
    hist[3] = Yhistory;
    hist[4] = Ahistory;
    hist[5] = Bhistory;
    hist[6] = Mhistory;
    return hist;
}
