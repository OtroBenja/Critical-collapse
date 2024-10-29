
#include <stdlib.h>
#include "constants.h"

double** iteration(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration,double epsilon){

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
    //double *Phi_ko = malloc(sizeof(double)*nR);
    //double *Pi_ko = malloc(sizeof(double)*nR);
    double *r2 = malloc(sizeof(double)*nR);
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
    double _rk[4] = {1.0,0.5,0.5,1.0};

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
        r2[ir] = r[ir]*r[ir];
    }
    if(METRIC == 0){
        for(int i=1;i<nR;i++){
            //a[i] = 1.0 + 0.0001*pow(r[i],3)*pow(E,-pow((r[i]-10.0)/10,2));
            //alpha[i] = pow(E,r[i]/100);
        }
    }

    for(int i=0;i<iterations;i++){

        if(SAVE_MODE == 0){
            //Save values of Phi, Pi, phi, a and alpha
            if(save_count == save_iteration){
                //printf("iteration %d\n",i);
                if(MASS){
                    mass = 0;
                    for(int ir=0;ir<nR;ir++){
                        mass += 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r2[ir]/(a[ir]*a[ir])*deltaR;;
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
                        mass += 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r2[ir]/(a[ir]*a[ir])*deltaR;
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

        //Advance Pi and Phi using RK4 
        for(int ir=0;ir<nR;ir++){
            l_n[ir] = 0.0;
            k_n[ir] = 0.0;
            j_sum[ir] = 0.0;
            k_sum[ir] = 0.0;
            l_sum[ir] = 0.0;
        }
       //calculate jn, kn and ln
        for(int n=0;n<4;n++){
            for(int ir=0;ir<nR;ir++){
                Phi_rk[ir] = Phi[ir]+_rk[n]*deltaT*k_n[ir];
                 Pi_rk[ir] =  Pi[ir]+_rk[n]*deltaT*l_n[ir];
                Beta[ir] = 2.0*PI*r[ir]*((Pi_rk[ir])*(Pi_rk[ir]) + (Phi_rk[ir])*(Phi_rk[ir]));
            }
            for(int ir=0;ir<nR-1;ir++){
                Beta1_2[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*(
                    ( Pi_rk[ir] + Pi_rk[ir+1])*( Pi_rk[ir] + Pi_rk[ir+1]) +
                    (Phi_rk[ir] +Phi_rk[ir+1])*(Phi_rk[ir] +Phi_rk[ir+1]));
            }
            //First iterate the metric for this Runge-Kutta step
            alpha[0] = 1;
            for(int ir=0;ir<nR-1;ir++){
                //calculate m1 and n1
                m1 = deltaR*    a[ir]*(Beta[ir]-0.5*(a[ir]*a[ir]-1)/r[ir]);
                n1 = deltaR*alpha[ir]*(Beta[ir]+0.5*(a[ir]*a[ir]-1)/r[ir]);
                //calculate m2 and n2
                m2 = deltaR*    (a[ir]+0.5*m1)*(Beta1_2[ir]-0.5*((a[ir]+0.5*m1)*(a[ir]+0.5*m1)-1)/(r[ir]+0.5*deltaR));
                n2 = deltaR*(alpha[ir]+0.5*n1)*(Beta1_2[ir]+0.5*((a[ir]+0.5*m1)*(a[ir]+0.5*m1)-1)/(r[ir]+0.5*deltaR));
                //calculate m3 and n3
                m3 = deltaR*    (a[ir]+0.5*m2)*(Beta1_2[ir]-0.5*((a[ir]+0.5*m2)*(a[ir]+0.5*m2)-1)/(r[ir]+0.5*deltaR));
                n3 = deltaR*(alpha[ir]+0.5*n2)*(Beta1_2[ir]+0.5*((a[ir]+0.5*m2)*(a[ir]+0.5*m2)-1)/(r[ir]+0.5*deltaR));
                //calculate m4 and n4
                m4 = deltaR*    (a[ir]+m3)*(Beta[ir+1]-0.5*((a[ir]+m3)*(a[ir]+m3)-1)/r[ir+1]);
                n4 = deltaR*(alpha[ir]+n3)*(Beta[ir+1]+0.5*((a[ir]+m3)*(a[ir]+m3)-1)/r[ir+1]);
                //Calculate next step for a and alpha
                a[ir+1]     = a[ir]     +(m1 +2.0*(m2+m3) +m4)/6.0;
                alpha[ir+1] = alpha[ir] +(n1 +2.0*(n2+n3) +n4)/6.0;
            }
            temp = 1.0/(alpha[nR-1]*a[nR-1]);
            for(int ir=0;ir<nR;ir++){
                alpha[ir]   = alpha[ir]*temp;
                Gamma[ir]   =        alpha[ir]*( Pi_rk[ir])/(a[ir]);
                Epsilon[ir] = r2[ir]*alpha[ir]*(Phi_rk[ir])/(a[ir]);
            }
            for(int ir=0;ir<5;ir++){
                rPhi[ir] = r[nR-5 +ir]*Phi_rk[nR-5 +ir];
                 rPi[ir] = r[nR-5 +ir]* Pi_rk[nR-5 +ir];
            }



            j_n[0] = Gamma[0];
            k_n[0] = 0;
            l_n[0] = alpha[0]*leftmost_D1(Phi_rk,0,deltaR)/a[0];
            j_n[1] = Gamma[1];
            k_n[1] = leftmid_D1(Gamma,1,deltaR);
            l_n[1] = leftmid_D1(Epsilon,1,deltaR)/r2[1];
            //#pragma omp parallel for
            for(int ir=2;ir<nR-2;ir++){
                j_n[ir] = Gamma[ir];
                k_n[ir] = centered_D1(Gamma,ir,deltaR);
                l_n[ir] = centered_D1(Epsilon,ir,deltaR)/r2[ir];
            }
            j_n[nR-2] = Gamma[nR-2];
            k_n[nR-2] = rightmid_D1(Gamma,nR-2,deltaR);
            l_n[nR-2] = rightmid_D1(Epsilon,nR-2,deltaR)/r2[nR-2];
            j_n[nR-1] = Gamma[nR-1];
            k_n[nR-1] = -rightmost_D1(rPhi,4,deltaR)/r[nR-1];
            l_n[nR-1] = -rightmost_D1( rPi,4,deltaR)/r[nR-1];

            //Calculate phi, Phi and Pi on next step
            for(int ir=0;ir<nR;ir++){
                j_sum[ir] += rk[n]*j_n[ir];
                k_sum[ir] += rk[n]*k_n[ir];
                l_sum[ir] += rk[n]*l_n[ir];
            }
        }

        //Calculate phi, Phi and Pi on next step
        for(int ir=0;ir<nR;ir++){
            phi[ir] += _6*deltaT*j_sum[ir];
            Phi[ir] += _6*deltaT*k_sum[ir];
            Pi[ir]  += _6*deltaT*l_sum[ir];
        }
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