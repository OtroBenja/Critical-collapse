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
//double **all_subgrids[10];
double norm_global;
double maxR_global;

#include "derivatives.h"
#include "initialize.h"
#include "integration.h"
#include "print.h"

//Function to integrate the whole metric
void full_metric(int level,int n_levels, double ***all_subgrids){
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
        //printf("grid %d; deltaR %lf; nR %d\n",i,deltaR,nR);

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

        //Iteraate metric and save final value to continue on the next grid
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
    double minR = *(grid_l[9]);
    double maxR = *(grid_l[10]);

    double* Beta = malloc(sizeof(double)*nR);
    double* Beta_p05 = malloc(sizeof(double)*nR);
    double* Gamma = malloc(sizeof(double)*nR);
    double* Epsilon = malloc(sizeof(double)*nR);
    double* Phi_rk = malloc(sizeof(double)*nR);
    double*  Pi_rk = malloc(sizeof(double)*nR);

    double* Jn = malloc(sizeof(double)*nR);
    double* Kn = malloc(sizeof(double)*nR);
    double* Ln = malloc(sizeof(double)*nR);
    double* Jsum = malloc(sizeof(double)*nR);
    double* Ksum = malloc(sizeof(double)*nR);
    double* Lsum = malloc(sizeof(double)*nR);

    //Get points from coarser grid (to the right)
    double   r_right[6];
    double phi_right[6];
    double Phi_right[6];
    double  Pi_right[6];
    double   a_right[6];
    double   alpha_right[6];
    double   Gamma_right[6];
    double Epsilon_right[6];
    double   *r_coarser = all_subgrids[grid_n-1][0];
    double *phi_coarser = all_subgrids[grid_n-1][1];
    double *Phi_coarser = all_subgrids[grid_n-1][2];
    double  *Pi_coarser = all_subgrids[grid_n-1][3];
    double   *a_coarser = all_subgrids[grid_n-1][4];
    double *alpha_coarser = all_subgrids[grid_n-1][5];
    for(int i=0;i<4;i++){
          r_right[i] =   r[nR-4+i];
        phi_right[i] = phi[nR-4+i];
        Phi_right[i] = Phi[nR-4+i];
         Pi_right[i] =  Pi[nR-4+i];
          a_right[i] =   a[nR-4+i];
        alpha_right[i] = alpha[nR-4+i];
    }
    //Copy or interpolate from coarser grid
      r_right[4] = 0.5*(  r_right[3] +   r_coarser[2]);
    phi_right[4] = 0.5*(phi_right[3] + phi_coarser[2]);
    Phi_right[4] = 0.5*(Phi_right[3] + Phi_coarser[2]);
     Pi_right[4] = 0.5*( Pi_right[3] +  Pi_coarser[2]);
      a_right[4] = 0.5*(  a_right[3] +   a_coarser[2]);
    alpha_right[4] = 0.5*(alpha_right[3] + alpha_coarser[2]);
      r_right[5] =   r_coarser[2];
    phi_right[5] = phi_coarser[2];
    Phi_right[5] = Phi_coarser[2];
     Pi_right[5] =  Pi_coarser[2];
      a_right[5] =   a_coarser[2];
    alpha_right[5] = alpha_coarser[2];

    //Get points from finer grid (to the left)
    double   r_left[6];
    double phi_left[6];
    double Phi_left[6];
    double  Pi_left[6];
    double   a_left[6];
    double   alpha_left[6];
    double   Gamma_left[6];
    double Epsilon_left[6];
    double   *r_finer = all_subgrids[grid_n+1][0];
    double *phi_finer = all_subgrids[grid_n+1][1];
    double *Phi_finer = all_subgrids[grid_n+1][2];
    double  *Pi_finer = all_subgrids[grid_n+1][3];
    double   *a_finer = all_subgrids[grid_n+1][4];
    double *alpha_finer = all_subgrids[grid_n+1][5];
    for(int i=0;i<4;i++){
          r_left[2+i] =   r[i];
        phi_left[2+i] = phi[i];
        Phi_left[2+i] = Phi[i];
         Pi_left[2+i] =  Pi[i];
          a_left[2+i] =   a[i];
        alpha_left[2+i] = alpha[i];
    }
    //Copy from finer grid
      r_left[0] =   r_finer[nR-1-2];
    phi_left[0] = phi_finer[nR-1-2];
    Phi_left[0] = Phi_finer[nR-1-2];
     Pi_left[0] =  Pi_finer[nR-1-2];
      a_left[0] =   a_finer[nR-1-2];
    alpha_left[0] = alpha_finer[nR-1-2];
      r_left[1] =   r_finer[nR-1-4];
    phi_left[1] = phi_finer[nR-1-4];
    Phi_left[1] = Phi_finer[nR-1-4];
     Pi_left[1] =  Pi_finer[nR-1-4];
      a_left[1] =   a_finer[nR-1-4];
    alpha_left[1] = alpha_finer[nR-1-4];

    //Calculate Gamma and Epsilon for both sides
    for(int ir=0;ir<6;ir++){
          Gamma_left[ir] =                       alpha_left[ir]*( Pi_left[ir])/(a_left[ir]);
        Epsilon_left[ir] = r_left[ir]*r_left[ir]*alpha_left[ir]*(Phi_left[ir])/(a_left[ir]);
          Gamma_right[ir] =                       alpha_right[ir]*( Pi_right[ir])/(a_right[ir]);
        Epsilon_right[ir] = r_right[ir]*r_right[ir]*alpha_right[ir]*(Phi_right[ir])/(a_right[ir]);
    }

    for(int it=0;it<nT;it++){
        //Advance Pi and Phi using RK4
        for(int ir=0;ir<nR;ir++){
              Ln[ir] = 0.0;
              Kn[ir] = 0.0;
            Jsum[ir] = 0.0;
            Ksum[ir] = 0.0;
            Lsum[ir] = 0.0;
        }
        for(int n=0;n<4;n++){
            //Calculate intermediate RK values for Phi, Pi and auxiliar variable Beta
            for(int ir=0;ir<nR;ir++){
                Phi_rk[ir] = Phi[ir]+_rk[n]*dT*Kn[ir];
                 Pi_rk[ir] =  Pi[ir]+_rk[n]*dT*Ln[ir];
                  Beta[ir] = 2.0*PI*r[ir]*((Pi_rk[ir])*(Pi_rk[ir]) + (Phi_rk[ir])*(Phi_rk[ir]));
            }
            for(int ir=0;ir<nR-1;ir++){
                Beta_p05[ir] = 0.5*PI*(r[ir]+0.5*dR)*(
                    ( Pi_rk[ir] + Pi_rk[ir+1])*( Pi_rk[ir] + Pi_rk[ir+1]) +
                    (Phi_rk[ir] +Phi_rk[ir+1])*(Phi_rk[ir] +Phi_rk[ir+1]));
            }
            //Iterate the metric for this Runge-Kutta step
            //metric_iteration(a0, alpha0, Beta, Beta_p05, a, alpha, r, nR, dR, norm_global);
            full_metric(1,N_LEVELS,all_subgrids);

            //Calculate auxiliar varaibles Gamma, Epsilon, rPhi and rPi
            for(int ir=0;ir<nR;ir++){
                Gamma[ir]   =             alpha[ir]*( Pi_rk[ir])/(a[ir]);
                Epsilon[ir] = r[ir]*r[ir]*alpha[ir]*(Phi_rk[ir])/(a[ir]);
            }
            for(int ir=0;ir<5;ir++){
                rPhi[ir] = r[nR-5 +ir]*Phi_rk[nR-5 +ir];
                 rPi[ir] = r[nR-5 +ir]* Pi_rk[nR-5 +ir];
            }

            //calculate jn, kn and ln
            Jn[0] = Gamma[0];
            Jn[1] = Gamma[1];
            if(minR==0){
                Kn[0] = 0;
                Ln[0] = alpha[0]*leftmost_D1(Phi_rk,0,dR)/a[0];
                Kn[1] = leftmid_D1(Gamma,1,dR);
                Ln[1] = leftmid_D1(Epsilon,1,dR)/(r[1]*r[1]);
            } else {
                Kn[0] = centered_D1(Gamma_left  ,2,dR);
                Ln[0] = centered_D1(Epsilon_left,2,dR)/(r_left[2]*r_left[2]);
                Kn[1] = centered_D1(Gamma_left  ,3,dR);
                Ln[1] = centered_D1(Epsilon_left,3,dR)/(r_left[3]*r_left[3]);
            }
            //#pragma omp parallel for
            for(int ir=2;ir<nR-2;ir++){
                Jn[ir] = Gamma[ir];
                Kn[ir] = centered_D1(Gamma,ir,dR);
                Ln[ir] = centered_D1(Epsilon,ir,dR)/(r[ir]*r[ir]);
            }
            Jn[nR-2] = Gamma[nR-2];
            Jn[nR-1] = Gamma[nR-1];

            if(maxR==maxR_global){
                Kn[nR-2] = rightmid_D1(Gamma,nR-2,dR);
                Ln[nR-2] = rightmid_D1(Epsilon,nR-2,dR)/(r[nR-2]*r[nR-2]);
                Kn[nR-1] = -rightmost_D1(rPhi,4,dR)/r[nR-1];
                Ln[nR-1] = -rightmost_D1( rPi,4,dR)/r[nR-1];
            } else {
                Kn[nR-2] = centered_D1(Gamma_right  ,2,dR);
                Ln[nR-2] = centered_D1(Epsilon_right,2,dR)/(r_right[2]*r_right[2]);
                Kn[nR-1] = centered_D1(Gamma_right  ,3,dR);
                Ln[nR-1] = centered_D1(Epsilon_right,3,dR)/(r_right[3]*r_right[3]);
            }

            for(int ir=0;ir<nR;ir++){
                Jsum[ir] += rk[n]*Jn[ir];
                Ksum[ir] += rk[n]*Kn[ir];
                Lsum[ir] += rk[n]*Ln[ir];
            }
        }
        
        //Calculate phi, Phi and Pi on next step
        for(int ir=0;ir<nR;ir++){
            phi[ir] += _6*dT*Jsum[ir];
            Phi[ir] += _6*dT*Ksum[ir];
             Pi[ir] += _6*dT*Lsum[ir];
        }
    }
    free(Beta);
    free(Beta_p05);
    free(Gamma);
    free(Epsilon);
    free(Phi_rk);
    free(Pi_rk);
    free(Jn);
    free(Kn);
    free(Ln);
    free(Jsum);
    free(Ksum);
    free(Lsum);

    if(grid_n>1){
        double **grid_next = all_subgrids[grid_n-1];
        double *phi_next = grid_next[1];
        double *Phi_next = grid_next[2];
        double  *Pi_next = grid_next[3];

        phi_next[0] = phi[nR-1];
        Phi_next[0] = Phi[nR-1];
         Pi_next[0] =  Pi[nR-1];
    }
}

double **initialize_subgrid(double fType,double *model_params,double initial_r, double final_r, double deltaR){
    double **initial_field;
    initial_field = initialize_field(fType,model_params,deltaR,initial_r,final_r);

    int nR = (int)((final_r-initial_r)/deltaR);
    double nR_f = (double)nR;
    double deltaT = deltaR/5.;
    double *a     = malloc(sizeof(double)*nR);
    double *alpha = malloc(sizeof(double)*nR);
    for(int ir=0;ir<nR;ir++){
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
    rf_p[0] = initial_field[0][nR-1];

    printf("nR: %d\n",nR);
    printf("dR_pointer: %lf\n",*dR_p);
    printf("dT_pointer: %lf\n",*dT_p);
    printf("grid range: %lf - %lf\n\n",*ri_p,*rf_p);

    double **subgrid = malloc(sizeof(double*)*11);
    subgrid[0]  = initial_field[0]; // r values of the grid
    subgrid[1]  = initial_field[1]; // phi values of the grid
    subgrid[2]  = initial_field[2]; // Phi values of the grid
    subgrid[3]  = initial_field[3]; // Pi values of the grid
    subgrid[4]  = a;     // a values of the grid
    subgrid[5]  = alpha; // alpha values of the grid
    subgrid[6]  = dR_p;  // deltaR of the grid (as pointer)
    subgrid[7]  = dT_p;  // deltaT of the grid (as pointer)
    subgrid[8]  = nR_p;  // number of points in the grid (as pointer)
    subgrid[9]  = ri_p;  // min r of the grid (as pointer)
    subgrid[10] = rf_p;  // max r of the grid (as pointer)

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
        else{minR_l = pow(2,-n)*maxR;}

        if(n==1){maxR_l = maxR+deltaR*1.01;}
        else{maxR_l = pow(2,-(n-1))*(maxR+deltaR*1.01);}
        double deltaR_l = pow(2,-(n-1))*deltaR;
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
    full_metric(1,n_levels,all_subgrids);

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
                    for(int ir_l=0; ir_l<nR_l-1; ir_l += pow(2,n_levels-(i+1))){
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
                    for(int ir_l=0; ir_l<nR_l-1; ir_l += pow(2,n_levels-(i+1))){
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
        full_metric(1,n_levels,all_subgrids);
        int count_level[n_levels-1];
        //bool trigger_level[n_levels];
        for(int i=0;i<n_levels;i++){
            count_level[i] = 0;
            //trigger_level[i] = false;
        }
        //trigger_level[0] = true;
        for(int i=0;i<n_levels-1;){
            //if(trigger_level[i]){
                integration_l(all_subgrids,n_levels-i,1,1,1);
                count_level[i]++;
                if(count_level[i]==2){
                    //trigger_level[i+1]=true;
                    count_level[i]=0;
                    i++;
                } else {
                    i = 0;
                }

            //}
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
    maxR_global = maxR;
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

