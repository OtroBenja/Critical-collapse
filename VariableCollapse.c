#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#if defined(_OPENMP)
#include <omp.h>
#else 
int omp_get_num_threads(){return 1;}
int omp_get_max_threads(){return 1;}
int omp_get_thread_num(){return 1;}
void omp_set_num_threads(){return;}
#endif

#include "derivatives.h"
#include "initialize.h"
#include "metric.h"
#include "print.h"


void integration_l(double *r_l, double *phi_l, double *Phi_l, double *Pi_l, double *a_l, double *alpha_l,
                    int nR_l, double dR_l, int nT_l, double dT_l){
    double  rk[4] = {1.0,2.0,2.0,1.0};
    double _rk[4] = {1.0,0.5,0.5,1.0};
    double *rPhi_l = malloc(sizeof(double)*5);
    double  *rPi_l = malloc(sizeof(double)*5);
    double temp1;
    double temp2;

    double* Beta_l = malloc(sizeof(double)*nR_l);
    double* Beta_p05_l = malloc(sizeof(double)*nR_l);
    double* Gamma_l = malloc(sizeof(double)*nR_l);
    double* Epsilon_l = malloc(sizeof(double)*nR_l);
    double* Phi_rk_l = malloc(sizeof(double)*nR_l);
    double*  Pi_rk_l = malloc(sizeof(double)*nR_l);

    double* Jn_l = malloc(sizeof(double)*nR_l);
    double* Kn_l = malloc(sizeof(double)*nR_l);
    double* Ln_l = malloc(sizeof(double)*nR_l);
    double* Jsum_l = malloc(sizeof(double)*nR_l);
    double* Ksum_l = malloc(sizeof(double)*nR_l);
    double* Lsum_l = malloc(sizeof(double)*nR_l);
    for(int it=0;it<nT_l;it++){
        //Advance Pi and Phi using RK4
        for(int ir=0;ir<nR_l;ir++){
            Ln_l[ir] = 0.0;
            Kn_l[ir] = 0.0;
            Jsum_l[ir] = 0.0;
            Ksum_l[ir] = 0.0;
            Lsum_l[ir] = 0.0;
        }
        //calculate jn, kn and ln
        for(int n=0;n<4;n++){
            for(int ir=0;ir<nR_l;ir++){
                Phi_rk_l[ir] = Phi_l[ir]+_rk[n]*dT_l*Kn_l[ir];
                 Pi_rk_l[ir] =  Pi_l[ir]+_rk[n]*dT_l*Ln_l[ir];
                Beta_l[ir] = 2.0*PI*r_l[ir]*((Pi_rk_l[ir])*(Pi_rk_l[ir]) + (Phi_rk_l[ir])*(Phi_rk_l[ir]));
            }
            for(int ir=0;ir<nR_l-1;ir++){
                Beta_p05_l[ir] = 0.5*PI*(r_l[ir]+0.5*dR_l)*(
                    ( Pi_rk_l[ir] + Pi_rk_l[ir+1])*( Pi_rk_l[ir] + Pi_rk_l[ir+1]) +
                    (Phi_rk_l[ir] +Phi_rk_l[ir+1])*(Phi_rk_l[ir] +Phi_rk_l[ir+1]));
            }
            //printf("print_l 3\n");
            //First iterate the metric for this Runge-Kutta step
            //printf("print_l 4\n");
            metric_iteration(Beta_l, Beta_p05_l, a_l, alpha_l, r_l, nR_l, dR_l);

            for(int ir=0;ir<nR_l;ir++){
                Gamma_l[ir]   =                 alpha_l[ir]*( Pi_rk_l[ir])/(a_l[ir]);
                Epsilon_l[ir] = r_l[ir]*r_l[ir]*alpha_l[ir]*(Phi_rk_l[ir])/(a_l[ir]);
            }
            for(int ir=0;ir<5;ir++){
                rPhi_l[ir] = r_l[nR_l-5 +ir]*Phi_rk_l[nR_l-5 +ir];
                 rPi_l[ir] = r_l[nR_l-5 +ir]* Pi_rk_l[nR_l-5 +ir];
            }

            Jn_l[0] = Gamma_l[0];
            Kn_l[0] = 0;
            Ln_l[0] = alpha_l[0]*leftmost_D1(Phi_rk_l,0,dR_l)/a_l[0];
            Jn_l[1] = Gamma_l[1];
            Kn_l[1] = leftmid_D1(Gamma_l,1,dR_l);
            Ln_l[1] = leftmid_D1(Epsilon_l,1,dR_l)/(r_l[1]*r_l[1]);
            //#pragma omp parallel for
            for(int ir=2;ir<nR_l-2;ir++){
                Jn_l[ir] = Gamma_l[ir];
                Kn_l[ir] = centered_D1(Gamma_l,ir,dR_l);
                Ln_l[ir] = centered_D1(Epsilon_l,ir,dR_l)/(r_l[ir]*r_l[ir]);
            }
            Jn_l[nR_l-2] = Gamma_l[nR_l-2];
            Kn_l[nR_l-2] = rightmid_D1(Gamma_l,nR_l-2,dR_l);
            Ln_l[nR_l-2] = rightmid_D1(Epsilon_l,nR_l-2,dR_l)/(r_l[nR_l-2]*r_l[nR_l-2]);
            Jn_l[nR_l-1] = Gamma_l[nR_l-1];
            Kn_l[nR_l-1] = -rightmost_D1(rPhi_l,4,dR_l)/r_l[nR_l-1];
            Ln_l[nR_l-1] = -rightmost_D1( rPi_l,4,dR_l)/r_l[nR_l-1];

            for(int ir=0;ir<nR_l;ir++){
                Jsum_l[ir] += rk[n]*Jn_l[ir];
                Ksum_l[ir] += rk[n]*Kn_l[ir];
                Lsum_l[ir] += rk[n]*Ln_l[ir];
            }
        }
        
        //Calculate phi, Phi and Pi on next step
        for(int ir=0;ir<nR_l;ir++){
            phi_l[ir] += _6*dT_l*Jsum_l[ir];
            Phi_l[ir] += _6*dT_l*Ksum_l[ir];
             Pi_l[ir] += _6*dT_l*Lsum_l[ir];
        }
    }
}

double **variable_iteration(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration,double epsilon){
    int nR = maxR/deltaR;
    int noSaveNR = MIN_R/deltaR;
    int saveNR = (maxR-MIN_R)/deltaR;
    double deltaT = deltaR/5.;
    double mass;
    double *a = malloc(sizeof(double)*nR);
    double *alpha = malloc(sizeof(double)*nR);
    double *Beta = malloc(sizeof(double)*nR);
    double *Beta1_2 = malloc(sizeof(double)*nR);
    double *Gamma = malloc(sizeof(double)*nR);
    double *Epsilon = malloc(sizeof(double)*nR);
    double *Phi_rk = malloc(sizeof(double)*nR);
    double * Pi_rk = malloc(sizeof(double)*nR);
    double m1, n1, m2, n2, m3, n3, m4, n4;

    double *rPhi = malloc(sizeof(double)*5);
    double *rPi = malloc(sizeof(double)*5);
    double *j_n = malloc(sizeof(double)*nR);
    double *k_n = malloc(sizeof(double)*nR);
    double *l_n = malloc(sizeof(double)*nR);
    double *j_sum = malloc(sizeof(double)*nR);
    double *k_sum = malloc(sizeof(double)*nR);
    double *l_sum = malloc(sizeof(double)*nR);
    double  rk[4] = {1.0,2.0,2.0,1.0};
    double _rk[4] = {0.0,0.5,0.5,1.0};

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

    for(int ir=0;ir<nR;ir++){
        a[ir] = 1.;
        alpha[ir] = 1.;
    }

    printf("print1\n");
    for(int i=0;i<iterations;i++){

        if(SAVE_MODE == 0){
            //Save values of Phi, Pi, phi, a and alpha
            if(save_count == save_iteration){
                //printf("iteration %d\n",i);
                if(MASS){
                    mass = 0;
                    for(int ir=0;ir<nR;ir++){
                        mass += 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir])*deltaR;
                    }
                    Mhistory[i/save_iteration] = mass;
                }
                for(int ir=0;ir<(nR/SAVE_RES);ir++){
                    //Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]*Phi[ir*SAVE_RES]*sqrt(2*PI)/a[ir*SAVE_RES];
                    //Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]* Pi[ir*SAVE_RES]*sqrt(2*PI)/a[ir*SAVE_RES];
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
                    mass = 0;
                    for(int ir=0;ir<nR;ir++){
                        mass += 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir])*deltaR;
                    }
                    Mhistory[save_iter] = mass;
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
        integration_l(r,phi,Phi,Pi,a,alpha,nR,deltaR,1,deltaT);
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

    initial_conditions = initialize_field(fType,model_parameters,deltaR,maxR);
    r = initial_conditions[0];
    phi = initial_conditions[1];
    Phi = initial_conditions[2];
    Pi = initial_conditions[3];

    //Pass initial conditions to iteration
    if((argc>9) && atoi(argv[9])) omp_set_num_threads(atoi(argv[9]));
    time_t initTime = time(NULL);
    hist = variable_iteration(r,phi,Phi,Pi,deltaR,maxR,iterations,SAVE_ITERATION,epsilon);
    time_t finalTime = time(NULL);
    int nP = omp_get_max_threads();
    time_t timeDelta = (finalTime-initTime);

    //Print simulation history to a file
    if(SAVE_MODE == 0)
        print_data(hist,fType,model_parameters,iterations,maxR,deltaR,nP,timeDelta,epsilon);
    else if(SAVE_MODE == 1){
        print_data(hist,fType,model_parameters,iterations-FIRST_ITERATION,maxR-MIN_R,deltaR,nP,timeDelta,epsilon);
        }
    printf("Finished, total time: %lds\n", timeDelta);
}




