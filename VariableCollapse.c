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
double **all_subgrids[10];

#include "derivatives.h"
#include "initialize.h"
#include "integration.h"
#include "print.h"


void full_metric(int level, int n_levels){
    double a0 = 1.0;
    double alpha0 = 1.0;
    double *a;
    double *alpha;
    for (int i=1;i<n_levels+1;i++){
        a = all_subgrids[i]
        = ;
    }
}

void integration_l(double **grid_l, double a0, double alpha0, int nT){
    double  rk[4] = {1.0,2.0,2.0,1.0};
    double _rk[4] = {1.0,0.5,0.5,1.0};
    double *rPhi = malloc(sizeof(double)*5);
    double  *rPi = malloc(sizeof(double)*5);
    double temp1;
    double temp2;

    double *r = grid_l[0];
    double *phi = grid_l[1];
    double *Phi = grid_l[2];
    double *Pi = grid_l[3];
    double *a = grid_l[4];
    double *alpha = grid_l[5];

    double dR = *(grid_l[6]);
    double dT = *(grid_l[7]);
    int nR = (int)(*(grid_l[8]));

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
    for(int it=0;it<nT;it++){
        //Advance Pi and Phi using RK4
        for(int ir=0;ir<nR;ir++){
            Ln[ir] = 0.0;
            Kn[ir] = 0.0;
            Jsum[ir] = 0.0;
            Ksum[ir] = 0.0;
            Lsum[ir] = 0.0;
        }
        //calculate jn, kn and ln
        for(int n=0;n<4;n++){
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
            //First iterate the metric for this Runge-Kutta step
            metric_iteration(a0, alpha0, Beta, Beta_p05, a, alpha, r, nR, dR, true);

            for(int ir=0;ir<nR;ir++){
                Gamma[ir]   =                 alpha[ir]*( Pi_rk[ir])/(a[ir]);
                Epsilon[ir] = r[ir]*r[ir]*alpha[ir]*(Phi_rk[ir])/(a[ir]);
            }
            for(int ir=0;ir<5;ir++){
                rPhi[ir] = r[nR-5 +ir]*Phi_rk[nR-5 +ir];
                 rPi[ir] = r[nR-5 +ir]* Pi_rk[nR-5 +ir];
            }

            Jn[0] = Gamma[0];
            Kn[0] = 0;
            Ln[0] = alpha[0]*leftmost_D1(Phi_rk,0,dR)/a[0];
            Jn[1] = Gamma[1];
            Kn[1] = leftmid_D1(Gamma,1,dR);
            Ln[1] = leftmid_D1(Epsilon,1,dR)/(r[1]*r[1]);
            //#pragma omp parallel for
            for(int ir=2;ir<nR-2;ir++){
                Jn[ir] = Gamma[ir];
                Kn[ir] = centered_D1(Gamma,ir,dR);
                Ln[ir] = centered_D1(Epsilon,ir,dR)/(r[ir]*r[ir]);
            }
            Jn[nR-2] = Gamma[nR-2];
            Kn[nR-2] = rightmid_D1(Gamma,nR-2,dR);
            Ln[nR-2] = rightmid_D1(Epsilon,nR-2,dR)/(r[nR-2]*r[nR-2]);
            Jn[nR-1] = Gamma[nR-1];
            Kn[nR-1] = -rightmost_D1(rPhi,4,dR)/r[nR-1];
            Ln[nR-1] = -rightmost_D1( rPi,4,dR)/r[nR-1];

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
}

double **initialize_subgrid(double fType,double *model_params,double initial_r, double final_r, double deltaR){
    double **initial_field;
    initial_field = initialize_field(fType,model_params,deltaR,initial_r,final_r);

    int nR = (int)((final_r-initial_r)/deltaR);
    double nR_f = (double)nR;
    double deltaT = deltaR/5.;
    double *a     = malloc(sizeof(double)*nR);
    double *alpha = malloc(sizeof(double)*nR);
    double *nR_p = malloc(sizeof(double));
    double *dR_p = malloc(sizeof(double));
    double *dT_p = malloc(sizeof(double));
    double *ri_p = malloc(sizeof(double));
    double *rf_p = malloc(sizeof(double));
    dR_p[0] = deltaR;
    dT_p[0] = deltaT;
    nR_p[0] = nR_f;
    ri_p[0] = initial_r;
    rf_p[0] = final_r;

    printf("nR: %d\n",nR);
    printf("nR_f: %lf\n",nR_f);
    printf("nR_pointer: %lf\n",*nR_p);
    printf("dR_pointer: %lf\n",*dR_p);
    printf("dT_pointer: %lf\n",*dT_p);

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

double **variable_iteration(int fType,double *model_params,double deltaR,int maxR,int iterations,int save_iteration){
    int noSaveNR = MIN_R/deltaR;
    int nR = (int)(maxR/deltaR);
    int saveNR = (maxR-MIN_R)/deltaR;
    double deltaT = deltaR/5.;
    double mass;

    double max_a = pow(TOLERANCE,-0.5); // Max value of 'a' according to g^rr tolerance
    bh_mass = 0;
    bh_radius = 0;
    is_bh = false;
    last_iteration = iterations;

    //Initialize main grid
    double **grid0;
    grid0 = initialize_subgrid(fType,model_params,0.0,maxR,deltaR);
    all_subgrids[0] = grid0;

    //Initialize all subgrids
    double **grid1;
    grid1 = initialize_subgrid(fType,model_params,0.5*maxR,maxR,deltaR);
    all_subgrids[1] = grid1;
    double **grid2;
    grid2 = initialize_subgrid(fType,model_params,0.25*maxR,0.5*(maxR+deltaR),0.5*deltaR);
    all_subgrids[2] = grid2;
    double **grid3;
    grid3 = initialize_subgrid(fType,model_params,0.125*maxR,0.25*(maxR+deltaR),0.25*deltaR);
    all_subgrids[3] = grid3;
    double **grid4;
    grid4 = initialize_subgrid(fType,model_params,0.0,0.125*(maxR+deltaR),0.125*deltaR);
    all_subgrids[4] = grid4;

    double temp;
    double *temp_pointer;
    double *Rhistory, *Fhistory, *Xhistory, *Yhistory, *Ahistory, *Bhistory, *Mhistory;
    if(SAVE_MODE == 0){
        Rhistory = malloc(sizeof(double)*(nR/SAVE_RES));
        Fhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
        Xhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
        Yhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
        Ahistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
        Bhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
        Mhistory = malloc(sizeof(double)*(iterations/save_iteration));
    }
    else if(SAVE_MODE == 1){
        Rhistory = malloc(sizeof(double)*saveNR);
        Fhistory = malloc(sizeof(double)*saveNR*(iterations-FIRST_ITERATION));
        Xhistory = malloc(sizeof(double)*saveNR*(iterations-FIRST_ITERATION));
        Yhistory = malloc(sizeof(double)*saveNR*(iterations-FIRST_ITERATION));
        Ahistory = malloc(sizeof(double)*saveNR*(iterations-FIRST_ITERATION));
        Bhistory = malloc(sizeof(double)*saveNR*(iterations-FIRST_ITERATION));
        Mhistory = malloc(sizeof(double)*(iterations-FIRST_ITERATION));
    }
    double** hist = malloc(sizeof(double*)*7);
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
    metric_iteration(1.0,1.0,Beta, Beta1_2, a, alpha, r, nR, deltaR, true);
    free(Beta);
    free(Beta1_2);
    printf("iteration stato\n");
    for(int i=0;i<iterations;i++){

        if(SAVE_MODE == 0){
            //Save values of Phi, Pi, phi, a and alpha
            if(save_count == save_iteration){
                //printf("iteration %d\n",i);
                if(MASS){
                    Mhistory[i/save_iteration] = get_mass(r,Phi,Pi,a,maxR,deltaR);;
                }
                for(int ir=0;ir<(nR/SAVE_RES);ir++){
                    Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =   Phi[ir*SAVE_RES];
                    Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =    Pi[ir*SAVE_RES];
                    Fhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =   phi[ir*SAVE_RES];
                    Ahistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =     a[ir*SAVE_RES];
                    Bhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = alpha[ir*SAVE_RES];
                }
                save_count=0;
            }
            save_count+=1;
        }

        if(SAVE_MODE == 1){
            //Save values of Phi, Pi, phi, a and alpha
            if(i >= FIRST_ITERATION){
                int save_iter = i - FIRST_ITERATION;
                //printf("iteration %d\n",i);
                if(MASS){
                    Mhistory[save_iter] = get_mass(r,Phi,Pi,a,maxR,deltaR);
                }
                for(int ir=0;ir<saveNR;ir++){
                    Xhistory[save_iter*saveNR + ir] =   Phi[noSaveNR + ir];
                    Yhistory[save_iter*saveNR + ir] =    Pi[noSaveNR + ir];
                    Fhistory[save_iter*saveNR + ir] =   phi[noSaveNR + ir];
                    Ahistory[save_iter*saveNR + ir] =     a[noSaveNR + ir];
                    Bhistory[save_iter*saveNR + ir] = alpha[noSaveNR + ir];
                }
            }
        }

        printf("print1\n");
        //Check if a collapse has happened if a > 1/TOLERANCE^2
        for(int ir=0;ir<nR;ir++){
            if(a[ir]>max_a){
                bh_radius = r[ir];
                bh_mass = get_mass(r,Phi,Pi,a,bh_radius,deltaR);
                last_iteration = i;
                is_bh = true;
                break;
            }
        }
        if(is_bh) break;
        printf("print2\n");

        //Integrate from smallest to bigger grid
        integration_l(grid0,1.0,1.0,1);
        printf("print3\n");
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
    double** hist;
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

    //Set Kreiss Oliger coefficient for dampening
    double epsilon = EPSILON;
    if((argc>5) && atof(argv[5])) epsilon  = atof(argv[5]);
    
    //Define simulation limits
    double deltaR = 0.01;
    if((argc>6) && atof(argv[6])) deltaR = atof(argv[6]);
    int maxR = 50;
    if((argc>7) && atoi(argv[7])) maxR = atoi(argv[7]);
    int iterations = ITERATIONS;
    if((argc>8) && atoi(argv[8])) iterations = atoi(argv[8]);
    printf("Total iterations: %d\n",iterations);

    //Pass initial conditions to iteration
    if((argc>9) && atoi(argv[9])) omp_set_num_threads(atoi(argv[9]));
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

