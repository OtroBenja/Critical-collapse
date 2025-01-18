#pragma once

#include <stdlib.h>
#include "constants.h"
#include "derivatives.h"

//Iterate on subgrid l, once
void subgrid_iteration(double ***all_subgrids,int grid_n, double a0, double alpha0, int count_level){
    double  rk[4] = {1.0,2.0,2.0,1.0};
    double _rk[4] = {1.0,0.5,0.5,1.0};
    double *rPhi = malloc(sizeof(double)*5);
    double  *rPi = malloc(sizeof(double)*5);
    double temp1;
    double temp2;
    //double ko_c = pow(-1,3)*EPSILON/5.;

    double **grid_l = all_subgrids[grid_n];
    double *r = grid_l[0];
    double *phi = grid_l[1];
    double *Phi = grid_l[2];
    double *Pi = grid_l[3];
    double *a = grid_l[4];
    double *alpha = grid_l[5];

    double dR = *(grid_l[6]);
    double dT = *(grid_l[7]);
    int nR = (int)(*(grid_l[8]));
    int nR_ex = nR + 2*GHOST_SIZE;
    double minR = *(grid_l[9]);
    double maxR = *(grid_l[10]);

    //double* Beta = malloc(sizeof(double)*nR_ex);
    //double* Beta_p05 = malloc(sizeof(double)*nR_ex);
    double* Gamma = malloc(sizeof(double)*nR_ex);
    double* Epsilon = malloc(sizeof(double)*nR_ex);
    //double *Phi_ko = malloc(sizeof(double)*nR);
    //double  *Pi_ko = malloc(sizeof(double)*nR);
    double* Phi_rk = malloc(sizeof(double)*nR_ex);
    double*  Pi_rk = malloc(sizeof(double)*nR_ex);
    //Beta     += GHOST_SIZE;
    //Beta_p05 += GHOST_SIZE;
    Gamma    += GHOST_SIZE;
    Epsilon  += GHOST_SIZE;
    Phi_rk   += GHOST_SIZE;
    Pi_rk    += GHOST_SIZE;

    double* Jn = malloc(sizeof(double)*nR_ex);
    double* Kn = malloc(sizeof(double)*nR_ex);
    double* Ln = malloc(sizeof(double)*nR_ex);
    double* Jsum = malloc(sizeof(double)*nR_ex);
    double* Ksum = malloc(sizeof(double)*nR_ex);
    double* Lsum = malloc(sizeof(double)*nR_ex);
    Jn   += GHOST_SIZE;
    Kn   += GHOST_SIZE;
    Ln   += GHOST_SIZE;
    Jsum += GHOST_SIZE;
    Ksum += GHOST_SIZE;
    Lsum += GHOST_SIZE;

    bool   load_finer = !((bool)(*(grid_l[11])));
    bool load_coarser = !((bool)(*(grid_l[12])));
    //if(minR==0.0) load_finer = false;
    //if(maxR==maxR_global) load_coarser = false;
    int left_idx = - GHOST_SIZE;
    int right_idx = nR + GHOST_SIZE;

    //load_ghostzones(all_subgrids,grid_l,grid_n);
    
    for(int ir=left_idx;ir<right_idx;ir++){
          Ln[ir] = 0.0;
          Kn[ir] = 0.0;
        Jsum[ir] = 0.0;
        Ksum[ir] = 0.0;
        Lsum[ir] = 0.0;
    }

    //Advance Pi and Phi using RK4
    for(int n=0;n<4;n++){

        //if(load_coarser){right_idx = nR + (GHOST_SIZE*(GRID_RATIO-count_level) - 2*n);
        if(load_coarser){right_idx = nR + (GHOST_SIZE - 2*n);
        } else {right_idx = nR;}
        //if(load_finer  ){left_idx  =    - (GHOST_SIZE*(GRID_RATIO-count_level) - 2*n);
        if(load_finer  ){left_idx  =    - (GHOST_SIZE - 2*n);
        } else {left_idx  = 0; }

        //Calculate intermediate RK values for Phi, Pi and auxiliar variable Beta
        for(int ir=left_idx; ir<right_idx; ir++){
            Phi_rk[ir] = Phi[ir]+_rk[n]*dT*Kn[ir];
             Pi_rk[ir] =  Pi[ir]+_rk[n]*dT*Ln[ir];
        }

        //Iterate the metric for this Runge-Kutta step
        //Calculate auxiliar variables Gamma, Epsilon, rPhi and rPi
        for(int ir=left_idx; ir<right_idx;ir++){
              Gamma[ir] =             alpha[ir]*( Pi_rk[ir])/(a[ir]);
            Epsilon[ir] = r[ir]*r[ir]*alpha[ir]*(Phi_rk[ir])/(a[ir]);
        }

        //calculate jn, kn and ln at non-boundaries
        for(int ir=left_idx;ir<right_idx;ir++)
            Jn[ir] = Gamma[ir];
        for(int ir=left_idx+2;ir<right_idx-2;ir++){
            Kn[ir] = centered_D1(Gamma,ir,dR);
            Ln[ir] = centered_D1(Epsilon,ir,dR)/(r[ir]*r[ir]);
        }

        //Calculate left boundary (origin) if neccesary
        if(!load_finer){
            Kn[0] = 0;
            Ln[0] = alpha[0]*leftmost_D1(Phi_rk,0,dR)/a[0];
            Kn[1] = leftmid_D1(Gamma,1,dR);
            Ln[1] = leftmid_D1(Epsilon,1,dR)/(r[1]*r[1]);
        }

        //Calculate right boundary (outgoing) if neccesary
        if(!load_coarser){
            for(int ir=0;ir<5;ir++){
                rPhi[ir] = r[nR-5 +ir]*Phi_rk[nR-5 +ir];
                 rPi[ir] = r[nR-5 +ir]* Pi_rk[nR-5 +ir];
            }
            Kn[nR-2] = rightmid_D1(Gamma,nR-2,dR);
            Ln[nR-2] = rightmid_D1(Epsilon,nR-2,dR)/(r[nR-2]*r[nR-2]);
            Kn[nR-1] = -rightmost_D1(rPhi,4,dR)/r[nR-1];
            Ln[nR-1] = -rightmost_D1( rPi,4,dR)/r[nR-1];
        }

        for(int ir=left_idx;ir<right_idx;ir++){
            Jsum[ir] += rk[n]*Jn[ir];
            Ksum[ir] += rk[n]*Kn[ir];
            Lsum[ir] += rk[n]*Ln[ir];
        }
    }
        
    //Calculate phi, Phi and Pi on next step
    for(int ir=left_idx;ir<right_idx;ir++){
        phi[ir] += _6*dT*Jsum[ir];
        Phi[ir] += _6*dT*Ksum[ir];
         Pi[ir] += _6*dT*Lsum[ir];
    }
    
    free(Gamma-GHOST_SIZE);
    free(Epsilon-GHOST_SIZE);
    free(Phi_rk-GHOST_SIZE);
    free(Pi_rk-GHOST_SIZE);
    free(Jn-GHOST_SIZE);
    free(Kn-GHOST_SIZE);
    free(Ln-GHOST_SIZE);
    free(Jsum-GHOST_SIZE);
    free(Ksum-GHOST_SIZE);
    free(Lsum-GHOST_SIZE);

    /*
    if(grid_n>1){
        double **grid_next = all_subgrids[grid_n-1];
        double *phi_next = grid_next[1];
        double *Phi_next = grid_next[2];
        double  *Pi_next = grid_next[3];

        phi_next[0] = phi[nR-1];
        Phi_next[0] = Phi[nR-1];
         Pi_next[0] =  Pi[nR-1];
    }
    */
}
