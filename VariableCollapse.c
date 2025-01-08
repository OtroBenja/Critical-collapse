#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#if defined(_OPENMP)
#include <omp.h>
#else 
int omp_get_num_threads(){return 1;}
int omp_get_max_threads(){return 1;}
int omp_get_thread_num(){return 1;}
void omp_set_num_threads(){return;}
#endif

double bh_mass;
double bh_radius;
bool is_bh;
int last_iteration;
double norm_global;
double maxR_global;

#include "derivatives.h"
#include "initialize.h"
#include "integration.h"
#include "print.h"

#define GHOST_SIZE 8

//Function to integrate the whole metric
void full_metric(int level,int n_levels, double ***all_subgrids){
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
        double     *Beta = malloc(sizeof(double)*nR);
        double *Beta_p05 = malloc(sizeof(double)*nR);
        for(int ir=0;ir<nR-1;ir++){
            Beta[ir] = 2.0*PI*r[ir]*((Pi[ir])*(Pi[ir]) + (Phi[ir])*(Phi[ir]));
            Beta_p05[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*(
                ( Pi[ir] + Pi[ir+1])*( Pi[ir] + Pi[ir+1]) +
                (Phi[ir] +Phi[ir+1])*(Phi[ir] +Phi[ir+1]));
        }
        Beta[nR-1] = 2.0*PI*r[nR-1]*((Pi[nR-1])*(Pi[nR-1]) + (Phi[nR-1])*(Phi[nR-1]));
        //Iterate metric and save final value to continue on the next grid
        //for(int j=0;j<10;j++){
            //printf("r%d: %e\n",j,r[j]);
            //printf("Phi%d: %lf; Pi%d %lf\n" ,j,Phi[j],j,Pi[j]);
            //printf("a%d: %lf; alpha%d %lf\n",j,  a[j],j,alpha[j]);
        //}
        metric_iteration(a0, alpha0, Beta, Beta_p05, a, alpha, r, nR, deltaR, 1.0);
        a0     = a[nR-1];
        alpha0 = alpha[nR-1];
        //printf("a0: %lf; alpha0 %lf\n\n",a0,alpha0);

        free(Beta);
        free(Beta_p05);
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

void copy_ghostzones(){}

//Iterate on subgrid l, nT times
void integration_l(double ***all_subgrids,int grid_n, double a0, double alpha0, int nT){
    double  rk[4] = {1.0,2.0,2.0,1.0};
    double _rk[4] = {1.0,0.5,0.5,1.0};
    double *rPhi = malloc(sizeof(double)*5);
    double  *rPi = malloc(sizeof(double)*5);
    double temp1;
    double temp2;

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

    double* Beta = malloc(sizeof(double)*nR_ex);
    double* Beta_p05 = malloc(sizeof(double)*nR_ex);
    double* Gamma = malloc(sizeof(double)*nR_ex);
    double* Epsilon = malloc(sizeof(double)*nR_ex);
    double* Phi_rk = malloc(sizeof(double)*nR_ex);
    double*  Pi_rk = malloc(sizeof(double)*nR_ex);
    Beta     += GHOST_SIZE;
    Beta_p05 += GHOST_SIZE;
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

    bool load_finer = true;
    bool load_coarser = true;
    if(minR==0.0) load_finer = false;
    if(maxR==maxR_global) load_coarser = false;

    int right_idx = nR + GHOST_SIZE;
    double   *Phi_coarse;
    double    *Pi_coarse;
    double     *a_coarse;
    double *alpha_coarse;
    if(load_coarser){
          Phi_coarse = all_subgrids[grid_n-1][2];
           Pi_coarse = all_subgrids[grid_n-1][3];
            a_coarse = all_subgrids[grid_n-1][4];
        alpha_coarse = all_subgrids[grid_n-1][5];
    }

    int left_idx = - GHOST_SIZE;
    double   *Phi_fine;
    double    *Pi_fine;
    double     *a_fine;
    double *alpha_fine;
    int nR_fine;
    if(load_finer){
          Phi_fine = all_subgrids[grid_n+1][2];
           Pi_fine = all_subgrids[grid_n+1][3];
            a_fine = all_subgrids[grid_n+1][4];
        alpha_fine = all_subgrids[grid_n+1][5];
           nR_fine = (int)(*(all_subgrids[grid_n+1][8]));
    }
    
    for(int it=0;it<nT;it++){

        for(int ir=left_idx;ir<right_idx;ir++){
              Ln[ir] = 0.0;
              Kn[ir] = 0.0;
            Jsum[ir] = 0.0;
            Ksum[ir] = 0.0;
            Lsum[ir] = 0.0;
        }
        if(load_coarser){
            //Get values from coarser grid (to the right)
            for(int ig=1;ig<=4;ig++){
                  Phi[nR-1 +GRID_RATIO*ig] =   Phi_coarse[ig];   
                   Pi[nR-1 +GRID_RATIO*ig] =    Pi_coarse[ig];   
                    a[nR-1 +GRID_RATIO*ig] =     a_coarse[ig];   
                alpha[nR-1 +GRID_RATIO*ig] = alpha_coarse[ig];   
            }
            for(int ig=0;ig<4;ig++){
            int idx = nR+GRID_RATIO*ig;
                  Phi[idx] = 0.5*(  Phi[idx-1]+  Phi[idx+1]);
                   Pi[idx] = 0.5*(   Pi[idx-1]+   Pi[idx+1]);
                    a[idx] = 0.5*(    a[idx-1]+    a[idx+1]);
                alpha[idx] = 0.5*(alpha[idx-1]+alpha[idx+1]);
            }
        }   
        if(load_finer){
            //Get values from finer grid (to the left)
            for(int ig=1;ig<=GHOST_SIZE;ig++){
                  Phi[-ig] =   Phi_fine[nR_fine-1 -GRID_RATIO*ig];
                   Pi[-ig] =    Pi_fine[nR_fine-1 -GRID_RATIO*ig];
                    a[-ig] =     a_fine[nR_fine-1 -GRID_RATIO*ig];
                alpha[-ig] = alpha_fine[nR_fine-1 -GRID_RATIO*ig];
            }
        }

        //Advance Pi and Phi using RK4
        for(int n=0;n<4;n++){

            if(load_coarser){right_idx = nR + (GHOST_SIZE - GRID_RATIO*n);
            } else {right_idx = nR;}
            if(load_finer  ){left_idx  =    - (GHOST_SIZE - GRID_RATIO*n);
            } else {left_idx  = 0; }

            //Calculate intermediate RK values for Phi, Pi and auxiliar variable Beta
            for(int ir=left_idx; ir<right_idx; ir++){
                Phi_rk[ir] = Phi[ir]+_rk[n]*dT*Kn[ir];
                 Pi_rk[ir] =  Pi[ir]+_rk[n]*dT*Ln[ir];
                  Beta[ir] = 2.0*PI*r[ir]*((Pi_rk[ir])*(Pi_rk[ir]) + (Phi_rk[ir])*(Phi_rk[ir]));
            }
            for(int ir=left_idx; ir<right_idx-1; ir++){
                Beta_p05[ir] = 0.5*PI*(r[ir]+0.5*dR)*(
                    ( Pi_rk[ir] + Pi_rk[ir+1])*( Pi_rk[ir] + Pi_rk[ir+1]) +
                    (Phi_rk[ir] +Phi_rk[ir+1])*(Phi_rk[ir] +Phi_rk[ir+1]));
            }

            //Iterate the metric for this Runge-Kutta step
            full_metric(1,N_LEVELS,all_subgrids);
            //Calculate auxiliar variables Gamma, Epsilon, rPhi and rPi
            for(int ir=left_idx; ir<right_idx;ir++){
                  Gamma[ir] =             alpha[ir]*( Pi_rk[ir])/(a[ir]);
                Epsilon[ir] = r[ir]*r[ir]*alpha[ir]*(Phi_rk[ir])/(a[ir]);
            }

            //calculate jn, kn and ln at non-boundaries
            for(int ir=left_idx+2;ir<right_idx-2;ir++){
                Jn[ir] = Gamma[ir];
                Kn[ir] = centered_D1(Gamma,ir,dR);
                Ln[ir] = centered_D1(Epsilon,ir,dR)/(r[ir]*r[ir]);
            }

            //Calculate left boundary (origin) if neccesary
            if(minR==0.0){
                Kn[0] = 0;
                Ln[0] = alpha[0]*leftmost_D1(Phi_rk,0,dR)/a[0];
                Kn[1] = leftmid_D1(Gamma,1,dR);
                Ln[1] = leftmid_D1(Epsilon,1,dR)/(r[1]*r[1]);
            }

            //Calculate right boundary (outgoing) if neccesary
            if(maxR==maxR_global){
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
    }
    free(Beta-GHOST_SIZE);
    free(Beta_p05-GHOST_SIZE);
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

double **initialize_subgrid(double fType,double *model_params,double initial_r, double final_r, double deltaR){
    double **initial_field;

    // If initial and final r are not in the boundaries of the simulation, add ghost zones

    double new_initial_r;
    int  left_ghost_size;
    if (initial_r == 0) {
        new_initial_r = initial_r;
        left_ghost_size = 0;
    } else {
        new_initial_r = initial_r - GHOST_SIZE*deltaR;
        left_ghost_size = GHOST_SIZE;
    }

    double   new_final_r;
    int right_ghost_size;
    if (final_r == maxR_global) {
        new_final_r = final_r + deltaR*0.01;
        right_ghost_size = 0;
    } else {
        new_final_r = final_r + GHOST_SIZE*deltaR + deltaR*0.01;
        right_ghost_size = GHOST_SIZE;
    }

    initial_field = initialize_field(fType,model_params,deltaR,new_initial_r,new_final_r);

    int nR = (int)((deltaR*0.01+final_r-initial_r)/deltaR);
    double nR_f = (double)nR;
    //int nR_ext = (int)((new_final_r-new_initial_r)/deltaR);
    int nR_ext = nR + left_ghost_size + right_ghost_size;
    double deltaT = deltaR/5.;
    double *a     = malloc(sizeof(double)*(nR_ext));
    double *alpha = malloc(sizeof(double)*(nR_ext));
    for(int ir=0; ir<(nR_ext); ir++){
        a[ir] = 1.0;
        alpha[ir] = 1.0;
    }
    double *nR_p = malloc(sizeof(double));
    double *dR_p = malloc(sizeof(double));
    double *dT_p = malloc(sizeof(double));
    double *ri_p = malloc(sizeof(double));
    double *rf_p = malloc(sizeof(double));
    dR_p[0] = deltaR;
    dT_p[0] = deltaT;
    nR_p[0] = nR_f;
    ri_p[0] = initial_r;
    //if(initial_r==0.0){
    //    initial_field[0][left_ghost_size] = 1.0E-50;
    //    ri_p[0] = 0.0;
    //} else {ri_p[0] = (initial_field[0]+left_ghost_size)[0];}
    //rf_p[0] = (initial_field[0]+left_ghost_size)[nR-1];
    rf_p[0] = final_r-deltaR;

    printf("nR: %d\n",nR);
    printf("nR_extended: %d\n",nR_ext);
    printf("dR_pointer: %lf\n",*dR_p);
    printf("dT_pointer: %lf\n",*dT_p);
    //printf("ri: %.20lf\n",initial_r);
    //printf("rf: %.20lf\n",final_r);
    printf("grid range: %lf - %lf\n\n",*ri_p,(initial_field[0]+left_ghost_size)[nR-1]);

    double **subgrid = malloc(sizeof(double*)*11);
    subgrid[0]  = initial_field[0]+left_ghost_size; // r values of the grid
    subgrid[1]  = initial_field[1]+left_ghost_size; // phi values of the grid
    subgrid[2]  = initial_field[2]+left_ghost_size; // Phi values of the grid
    subgrid[3]  = initial_field[3]+left_ghost_size; // Pi values of the grid
    subgrid[4]  =     a+left_ghost_size; // a values of the grid
    subgrid[5]  = alpha+left_ghost_size; // alpha values of the grid
    subgrid[6]  = dR_p; // deltaR of the grid (as pointer)
    subgrid[7]  = dT_p; // deltaT of the grid (as pointer)
    subgrid[8]  = nR_p; // number of points in the grid (as pointer)
    subgrid[9]  = ri_p; // min r of the grid (as pointer)
    subgrid[10] = rf_p; // max r of the grid (as pointer)

    return subgrid;
}

float **variable_iteration(int fType,double *model_params,double deltaR,double maxR,int iterations,int save_iteration){
    int noSaveNR = MIN_R/deltaR;
    int nR = (int)(maxR/deltaR);
    int saveNR = (maxR-MIN_R)/deltaR;
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
    grid0 = initialize_subgrid(fType,model_params,0.0,maxR+deltaR,deltaR);
    all_subgrids[0] = grid0;

    //Initialize all subgrids
    int n_levels = N_LEVELS;
    for(int n=1; n<n_levels+1; n++){
        double **gridN;
        double minR_l;
        double maxR_l;
        
        if(n==n_levels){minR_l = 0;}
        else{minR_l = pow(GRID_RATIO,-n)*maxR;}

        //if(n==1){maxR_l = maxR+deltaR*1.01;}
        if(n==1){maxR_l = maxR+deltaR;}
        else{maxR_l = pow(GRID_RATIO,-(n-1))*(maxR+deltaR);}
        //else{maxR_l = pow(GRID_RATIO,-(n-1))*(maxR+deltaR*1.01);}
        double deltaR_l = pow(GRID_RATIO,-(n-1))*deltaR;
        gridN = initialize_subgrid(fType,model_params,minR_l,maxR_l,deltaR_l);
        all_subgrids[n] = gridN;
    }

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
    full_metric(1,N_LEVELS,all_subgrids);

    printf("iteration stato\n");
    for(int i=0;i<iterations;i++){

        if(SAVE_MODE == 0){
            //Save values of Phi, Pi, phi, a and alpha
            if(save_count == save_iteration){
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
                    for(int ir_l=0; ir_l<nR_l-1; ir_l += pow(GRID_RATIO,n_levels-(i+1))){
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
                    for(int ir_l=0; ir_l<nR_l-1; ir_l += pow(GRID_RATIO,n_levels-(i+1))){
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

        //printf("print1\n");
        //Check if a collapse has happened if a > 1/TOLERANCE^2
        for(int ir=0;ir<nR;ir++){
            if(a[ir]>max_a){
                bh_radius = r[ir];
                mass = get_mass(r,Phi,Pi,a,maxR,deltaR);
                bh_mass = mass[ir];
                last_iteration = i;
                is_bh = true;
                break;
            }
        }
        if(is_bh) break;
        //printf("print2\n");
        //integration_l(grid0,1.0,1.0,1);
        //Integrate from smallest to bigger grid
        full_metric(1,N_LEVELS,all_subgrids);
        int count_level[n_levels-1];
        for(int i=0;i<n_levels;i++){
            count_level[i] = 0;
        }
        //trigger_level[0] = true;
        for(int i=0;i<n_levels-1;){
            integration_l(all_subgrids,n_levels-i,1,1,1);
            count_level[i]++;
            if(count_level[i]==GRID_RATIO){
                count_level[i]=0;
                i++;
            } else {
                i = 0;
            }
        }
        integration_l(all_subgrids,1,1,1,1);
        //printf("print3\n");
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


int main(int argc, char* argv[]){
    double* model_parameters = malloc(sizeof(double)*3);
    double** initial_conditions;
    float** hist;
    double* r;
    double* phi;
    double* Phi;
    double* Pi;
    //Define simulation parameters
    int fType = 0;
    if((argc>1) && atoi(argv[1])) fType = atoi(argv[1]);
    double p0 = 0.000008;
    if((argc>2) && atof(argv[2])) p0 = atof(argv[2]);
    double r0 = 20.;
    if((argc>3) && atof(argv[3])) r0 = atof(argv[3]);
    double d = 3.;
    if((argc>4) && atof(argv[4])) d  = atof(argv[4]);
    model_parameters[0] = p0;
    model_parameters[1] = r0;
    model_parameters[2] = d;
    
    //Define simulation limits
    double deltaR = 0.01;
    if((argc>5) && atof(argv[5])) deltaR = atof(argv[5]);
    double maxR = 50;
    if((argc>6) && atoi(argv[6])) maxR = atof(argv[6]);
    maxR_global = maxR+deltaR;
    int iterations = ITERATIONS;
    if((argc>7) && atoi(argv[7])) iterations = atoi(argv[7]);
    printf("Total iterations: %d\n",iterations);

    //Pass initial conditions to iteration
    if((argc>8) && atoi(argv[8])) omp_set_num_threads(atoi(argv[8]));
    time_t initTime = time(NULL);
    hist = variable_iteration(fType, model_parameters,deltaR,maxR,iterations,SAVE_ITERATION);
    time_t finalTime = time(NULL);
    int nP = omp_get_max_threads();
    time_t timeDelta = (finalTime-initTime);

    //Print simulation history to a file
    if(SAVE_MODE == 0)
        print_data(hist,fType,model_parameters,last_iteration,maxR,deltaR,nP,timeDelta);
    else if(SAVE_MODE == 1){
        print_data(hist,fType,model_parameters,last_iteration-FIRST_ITERATION,maxR-MIN_R,deltaR,nP,timeDelta);
        }
    printf("Finished, total time: %lds\n", timeDelta);
}

