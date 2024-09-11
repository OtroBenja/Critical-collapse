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

#define MASS 1
#define METRIC 2 // 0 = minkowski; 1 = choptuik; 2 = modified choptuik
#define SAVE_MODE 1 // 0 = save uniformly on every SAVE_RES and SAVE_ITERATION ; 1 = save all points after FIRST_ITERATION and with r > MIN_R
#define SAVE_RES 500
#define SAVE_ITERATION 100
#define FIRST_ITERATION 77500
#define MIN_R 45
#define ITERATIONS 80000
#define PI 3.141592653
#define E  2.718281828


double** initialize_field(int fType,double* model_parameters,double deltaR,int maxR){
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
    double q  = model_parameters[3];
    int nR = maxR/deltaR;
    double* r = malloc(sizeof(double)*nR);
    double* phi = malloc(sizeof(double)*nR);
    double* Phi = malloc(sizeof(double)*nR);
    double* Pi = malloc(sizeof(double)*nR);
    double** return_values = malloc(sizeof(double*)*4);

    //Define r and calculate initial phi
    r[0] = deltaR*0.1;
    for(int i=1;i<nR;i++)
        r[i] = i*deltaR;
    if(fType==0){
        for(int i=0;i<nR;i++)
            phi[i] = p0*tanh((r[i]-r0)/d);
    }
    if(fType==1){
        for(int i=0;i<nR;i++)
            phi[i] = p0*pow(r[i],3)*pow(E,-pow((r[i]-r0)/d,q));
    }
    if(fType==2){
        for(int i=0;i<nR;i++)
            phi[i] = p0;
    }
    if(fType==3){
        for(int i=0;i<nR;i++)
            phi[i] = p0*r[i];
    }

    //Calculate initial Phi
    Phi[0] = 0;
    Phi[1] = (phi[4] -6.0*phi[3] +18.0*phi[2]
        -10.0*phi[1] -3.0*phi[0])/(12.0*deltaR);
    for(int i=2;i<nR-2;i++){
        Phi[i] = (-phi[i+2] +8.0*phi[i+1] -8.0*phi[i-1] +phi[i-2])/(12.0*deltaR);
    }
    Phi[nR-2] = (-phi[nR-5] +6.0*phi[nR-4] -18.0*phi[nR-3] 
            +10.0*phi[nR-2] +3.0*phi[nR-1])/(12.0*deltaR);
    Phi[nR-1] = 0;

    //Set initial Pi to zero
    for(int i=0;i<nR;i++){
        Pi[i] = 0;
    }
        
    return_values[0] = r;
    return_values[1] = phi;
    return_values[2] = Phi;
    return_values[3] = Pi;
 
    return return_values;
}

double** iteration(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration){

    int nR = maxR/deltaR;
    int noSaveNR = MIN_R/deltaR;
    int saveNR = (maxR-MIN_R)/deltaR;
    double deltaT = deltaR/5.;
    double mass;
    double* a = malloc(sizeof(double)*nR);
    double* alpha = malloc(sizeof(double)*nR);
    double* Beta = malloc(sizeof(double)*nR);
    double* Beta1_2 = malloc(sizeof(double)*nR);
    double* Gamma = malloc(sizeof(double)*nR);
    double* Epsilon = malloc(sizeof(double)*nR);
    double* Zeta = malloc(sizeof(double)*nR);
    double* r2 = malloc(sizeof(double)*nR);
    double* j1 = malloc(sizeof(double)*nR);
    double* k1 = malloc(sizeof(double)*nR);
    double* l1 = malloc(sizeof(double)*nR);
    double* j2 = malloc(sizeof(double)*nR);
    double* k2 = malloc(sizeof(double)*nR);
    double* l2 = malloc(sizeof(double)*nR);
    double* j3 = malloc(sizeof(double)*nR);
    double* k3 = malloc(sizeof(double)*nR);
    double* l3 = malloc(sizeof(double)*nR);
    double* j4 = malloc(sizeof(double)*nR);
    double* k4 = malloc(sizeof(double)*nR);
    double* l4 = malloc(sizeof(double)*nR);
    double m1, n1, m2, n2, m3, n3, m4, n4;
    double temp;
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
        Gamma[ir] = 1.;
        r2[ir] = r[ir]*r[ir];
        Epsilon[ir] = r2[ir];
    }
    printf("Pre-iteration");
    for(int i=0;i<iterations;i++){
        
        //If on minkowski metric, do not solve a and alpha

        //Useful auxiliar variable
        #pragma omp parallel for
        for(int ir=0;ir<nR;ir++)
            Beta[ir] = 2.0*PI*r[ir]*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
        //If on choptuik metric, solve a and alpha with known values of PHI and PI using RK4
        if(METRIC == 1){
            #pragma omp parallel for
            for(int ir=0;ir<nR-1;ir++){
                Beta1_2[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*((Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])+(Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1]));
                //Beta1_2[ir] = 2.0*PI*(r[ir]+deltaR/2)*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
            }
            
            for(int ir=0;ir<nR-1;ir++){
                //calculate m1 and n1
                m1 =     a[ir]*(Beta[ir]-0.5*(a[ir]*a[ir]-1)/r[ir]);
                n1 = alpha[ir]*(Beta[ir]+0.5*(a[ir]*a[ir]-1)/r[ir]);
                //calculate m2 and n2
                m2 =     (a[ir]+0.5*deltaR*m1)*(Beta1_2[ir]-0.5*((a[ir]+0.5*deltaR*m1)*(a[ir]+0.5*deltaR*m1)-1)/(r[ir]+0.5*deltaR));
                n2 = (alpha[ir]+0.5*deltaR*n1)*(Beta1_2[ir]+0.5*((a[ir]+0.5*deltaR*m1)*(a[ir]+0.5*deltaR*m1)-1)/(r[ir]+0.5*deltaR));
                //calculate m3 and n3
                m3 =     (a[ir]+0.5*deltaR*m2)*(Beta1_2[ir]-0.5*((a[ir]+0.5*deltaR*m2)*(a[ir]+0.5*deltaR*m2)-1)/(r[ir]+0.5*deltaR));
                n3 = (alpha[ir]+0.5*deltaR*n2)*(Beta1_2[ir]+0.5*((a[ir]+0.5*deltaR*m2)*(a[ir]+0.5*deltaR*m2)-1)/(r[ir]+0.5*deltaR));
                //calculate m4 and n4
                m4 =     (a[ir]+deltaR*m3)*(Beta[ir+1]-0.5*((a[ir]+deltaR*m3)*(a[ir]+deltaR*m3)-1)/r[ir+1]);
                n4 = (alpha[ir]+deltaR*n3)*(Beta[ir+1]+0.5*((a[ir]+deltaR*m3)*(a[ir]+deltaR*m3)-1)/r[ir+1]);
                //Calculate next step for a and alpha
                a[ir+1]     = a[ir]     + deltaR*(m1 +2.0*(m2+m3) +m4)/6.0;
                alpha[ir+1] = alpha[ir] + deltaR*(n1 +2.0*(n2+n3) +n4)/6.0;
            }
            //Define auxiliar variables
            #pragma omp parallel for
            for(int ir=0;ir<nR;ir++){
                Gamma[ir] = alpha[ir]/a[ir];
                Epsilon[ir] = r2[ir]*Gamma[ir];
                Zeta[ir] = r[ir]/(a[ir]*a[ir]);
            }
        }
        //If on modified choptuik metric, solve s and sigma with known values of PHI and PI using RK4
        if(METRIC == 2){
            #pragma omp parallel for
            for(int ir=0;ir<nR-1;ir++){
                Beta1_2[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*((Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])+(Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1]));
                //Beta1_2[ir] = 2.0*PI*(r[ir]+deltaR/2)*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
            }
            
            for(int ir=0;ir<nR-1;ir++){
                //calculate m1 and n1
                temp = a[ir]*a[ir];
                m1 =     -a[ir]*(Beta[ir]-0.5*(1-temp)/(r[ir]*temp));
                n1 = -alpha[ir]*(Beta[ir]+0.5*(1-temp)/(r[ir]*temp));
                //calculate m2 and n2
                temp = (a[ir]+0.5*deltaR*m1)*(a[ir]+0.5*deltaR*m1);
                m2 =     -(a[ir]+0.5*deltaR*m1)*(Beta1_2[ir]-0.5*(1-temp)/((r[ir]+0.5*deltaR)*temp));
                n2 = -(alpha[ir]+0.5*deltaR*n1)*(Beta1_2[ir]+0.5*(1-temp)/((r[ir]+0.5*deltaR)*temp));
                //calculate m3 and n3
                temp = (a[ir]+0.5*deltaR*m2)*(a[ir]+0.5*deltaR*m2);
                m3 =     -(a[ir]+0.5*deltaR*m2)*(Beta1_2[ir]-0.5*(1-temp)/((r[ir]+0.5*deltaR)*temp));
                n3 = -(alpha[ir]+0.5*deltaR*n2)*(Beta1_2[ir]+0.5*(1-temp)/((r[ir]+0.5*deltaR)*temp));
                //calculate m4 and n4
                temp = (a[ir]+deltaR*m3)*(a[ir]+deltaR*m3);
                m4 =     -(a[ir]+deltaR*m3)*(Beta[ir+1]-0.5*(1-temp)/(r[ir+1]*temp));
                n4 = -(alpha[ir]+deltaR*n3)*(Beta[ir+1]+0.5*(1-temp)/(r[ir+1]*temp));
                //Calculate next step for a and alpha
                a[ir+1]     = a[ir]     + deltaR*(m1 +2.0*(m2+m3) +m4)/6.0;
                alpha[ir+1] = alpha[ir] + deltaR*(n1 +2.0*(n2+n3) +n4)/6.0;
            }
            //Define auxiliar variables
            #pragma omp parallel for
            for(int ir=0;ir<nR;ir++){
                Gamma[ir] = a[ir]/alpha[ir];
                Epsilon[ir] = r2[ir]*Gamma[ir];
                Zeta[ir] = r[ir]*a[ir]*a[ir];
            }
        }
        //If on choptuik metric 1st order, solve a and alpha with known values of PHI and PI using Euler's method
        if(METRIC == 3){
            
            for(int ir=0;ir<nR-1;ir++){
                //Calculate next step for a and alpha
                a[ir+1]     = a[ir]     + deltaR*    a[ir]*(Beta[ir]-0.5*(a[ir]*a[ir]-1)/r[ir]);
                alpha[ir+1] = alpha[ir] + deltaR*alpha[ir]*(Beta[ir]+0.5*(a[ir]*a[ir]-1)/r[ir]);
            }
            //Define auxiliar variables
            #pragma omp parallel for
            for(int ir=0;ir<nR;ir++){
                Gamma[ir] = alpha[ir]/a[ir];
                Epsilon[ir] = r2[ir]*Gamma[ir];
                Zeta[ir] = r[ir]/(a[ir]*a[ir]);
            }
        }
        //If on modified choptuik metric 1st order, solve s and sigma with known values of PHI and PI using Euler's method
        if(METRIC == 4){
            
            for(int ir=0;ir<nR-1;ir++){
                temp = a[ir]*a[ir];
                //Calculate next step for a and alpha
                a[ir+1]     = a[ir]     - deltaR*    a[ir]*(Beta[ir]-0.5*(1-temp)/(r[ir]*temp));
                alpha[ir+1] = alpha[ir] - deltaR*alpha[ir]*(Beta[ir]+0.5*(1-temp)/(r[ir]*temp));
            }
            //Define auxiliar variables
            #pragma omp parallel for
            for(int ir=0;ir<nR;ir++){
                Gamma[ir] = a[ir]/alpha[ir];
                Epsilon[ir] = r2[ir]*Gamma[ir];
                Zeta[ir] = r[ir]*a[ir]*a[ir];
            }
        }

        if(SAVE_MODE == 0){
            //Save values of Phi, Pi, phi, a and alpha
            if(save_count == save_iteration){
                //printf("iteration %d\n",i);
                if(MASS){
                    mass = 0;
                    for(int ir=0;ir<nR;ir++){
                        mass += Beta[ir]*Zeta[ir]*deltaR;
                    }
                    Mhistory[i/save_iteration] = mass;
                }
                #pragma omp parallel for
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
                        mass += Beta[ir]*Zeta[ir]*deltaR;
                    }
                    Mhistory[save_iter] = mass;
                }
                #pragma omp parallel for
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

        //calculate j1, k1 and l1
        j1[0] =          Gamma[0]* Pi[0];
        k1[0] =   (-25.0*Gamma[0]* Pi[0]   +48.0*Gamma[1]* Pi[1]   -36.0*Gamma[2]* Pi[2] 
                   +16.0*Gamma[3]* Pi[3]    -3.0*Gamma[4]* Pi[4])/(12.0*deltaR);
        l1[0] = (-25.0*Epsilon[0]*Phi[0] +48.0*Epsilon[1]*Phi[1] -36.0*Epsilon[2]*Phi[2]
                 +16.0*Epsilon[3]*Phi[3]  -3.0*Epsilon[4]*Phi[4])/(12.0*deltaR*r2[0]);
        j1[1] =    Gamma[1]* Pi[1];
        k1[1] =   (Gamma[4]* Pi[4]   -6.0*Gamma[3]* Pi[3]   +18.0*Gamma[2]* Pi[2]
             -10.0*Gamma[1]* Pi[1]   -3.0*Gamma[0]* Pi[0])/(12.0*deltaR);
        l1[1] = (Epsilon[4]*Phi[4] -6.0*Epsilon[3]*Phi[3] +18.0*Epsilon[2]*Phi[2]
           -10.0*Epsilon[1]*Phi[1] -3.0*Epsilon[0]*Phi[0])/(12.0*deltaR*r2[1]);
        #pragma omp parallel for
        for(int ir=2;ir<nR-2;ir++){
            j1[ir] =     Gamma[ir  ]* Pi[ir  ];
            k1[ir] =   (-Gamma[ir+2]* Pi[ir+2]   +8.0*Gamma[ir+1]* Pi[ir+1] 
                    -8.0*Gamma[ir-1]* Pi[ir-1]       +Gamma[ir-2]* Pi[ir-2])/(12.0*deltaR);
            l1[ir] = (-Epsilon[ir+2]*Phi[ir+2] +8.0*Epsilon[ir+1]*Phi[ir+1] 
                  -8.0*Epsilon[ir-1]*Phi[ir-1]     +Epsilon[ir-2]*Phi[ir-2])/(12.0*deltaR*r2[ir]);
        }
        j1[nR-2] =     Gamma[nR-2]* Pi[nR-2];
        k1[nR-2] = 0;
        //k1[nR-2] =   (-Gamma[nR-5]* Pi[nR-5]   +6.0*Gamma[nR-4]* Pi[nR-4]   -18.0*Gamma[nR-3]* Pi[nR-3] 
        //         +10.0*Gamma[nR-2]* Pi[nR-2]   +3.0*Gamma[nR-1]* Pi[nR-1])/(12.0*deltaR);
        k1[nR-2] = 0;
        //l1[nR-2] = (-Epsilon[nR-5]*Phi[nR-5] +6.0*Epsilon[nR-4]*Phi[nR-4] -18.0*Epsilon[nR-3]*Phi[nR-3] 
        //       +10.0*Epsilon[nR-2]*Phi[nR-2] +3.0*Epsilon[nR-1]*Phi[nR-1])/(12.0*deltaR*r2[nR-2]);
        j1[nR-1] =     Gamma[nR-1]* Pi[nR-1];
        k1[nR-1] = 0;
        //k1[nR-1] = phi[nR-1]/r2[nR-1] -Phi[nR-1]/r[nR-1]-0.5*(3*Phi[nR-1] -4*Phi[nR-2] +Phi[nR-3]);
        l1[nR-1] = 0;
        //l1[nR-1] = -2*Gamma[nR-1]*(phi[nR-1]+Gamma[nR-1]*Pi[nR-1]+0.5*r[nR-1]*Phi[nR-1])/r[nR-1] +
        //           -0.5*(3*Gamma[nR-1]*Gamma[nR-1]*Pi[nR-1] -4*Gamma[nR-2]*Gamma[nR-2]*Pi[nR-2] +Gamma[nR-3]*Gamma[nR-3]*Pi[nR-3]) +
        //           -0.5*phi[nR-1]*(3*Gamma[nR-1] -4*Gamma[nR-2] +Gamma[nR-3]) ;


        //calculate j2, k2 and l2
        j2[0] =          Gamma[0]*( Pi[0]+0.5*l1[0]*deltaT);
        k2[0] =   (-25.0*Gamma[0]*( Pi[0]+0.5*l1[0]*deltaT)   +48.0*Gamma[1]*( Pi[1]+0.5*l1[1]*deltaT)   -36.0*Gamma[2]*( Pi[2]+0.5*l1[2]*deltaT) 
                   +16.0*Gamma[3]*( Pi[3]+0.5*l1[3]*deltaT)    -3.0*Gamma[4]*( Pi[4]+0.5*l1[4]*deltaT))/(12.0*deltaR);
        l2[0] = (-25.0*Epsilon[0]*(Phi[0]+0.5*k1[0]*deltaT) +48.0*Epsilon[1]*(Phi[1]+0.5*k1[1]*deltaT) -36.0*Epsilon[2]*(Phi[2]+0.5*k1[2]*deltaT)
                 +16.0*Epsilon[3]*(Phi[3]+0.5*k1[3]*deltaT)  -3.0*Epsilon[4]*(Phi[4]+0.5*k1[4]*deltaT))/(12.0*deltaR*r2[0]);
        j2[1] =    Gamma[1]*( Pi[1]+0.5*l1[1]*deltaT);
        k2[1] =   (Gamma[4]*( Pi[4]+0.5*l1[4]*deltaT)   -6.0*Gamma[3]*( Pi[3]+0.5*l1[3]*deltaT)   +18.0*Gamma[2]*( Pi[2]+0.5*l1[2]*deltaT)
             -10.0*Gamma[1]*( Pi[1]+0.5*l1[1]*deltaT)   -3.0*Gamma[0]*( Pi[0]+0.5*l1[0]*deltaT))/(12.0*deltaR);
        l2[1] = (Epsilon[4]*(Phi[4]+0.5*k1[4]*deltaT) -6.0*Epsilon[3]*(Phi[3]+0.5*k1[3]*deltaT) +18.0*Epsilon[2]*(Phi[2]+0.5*k1[2]*deltaT)
           -10.0*Epsilon[1]*(Phi[1]+0.5*k1[1]*deltaT) -3.0*Epsilon[0]*(Phi[0]+0.5*k1[0]*deltaT))/(12.0*deltaR);
        #pragma omp parallel for
        for(int ir=2;ir<nR-2;ir++){
            j2[ir] =     Gamma[ir  ]*( Pi[ir  ]+0.5*l1[ir  ]*deltaT);
            k2[ir] =   (-Gamma[ir+2]*( Pi[ir+2]+0.5*l1[ir+2]*deltaT)   +8.0*Gamma[ir+1]*( Pi[ir+1]+0.5*l1[ir+1]*deltaT) 
                    -8.0*Gamma[ir-1]*( Pi[ir-1]+0.5*l1[ir-1]*deltaT)       +Gamma[ir-2]*( Pi[ir-2]+0.5*l1[ir-2]*deltaT))/(12.0*deltaR);
            l2[ir] = (-Epsilon[ir+2]*(Phi[ir+2]+0.5*k1[ir+2]*deltaT) +8.0*Epsilon[ir+1]*(Phi[ir+1]+0.5*k1[ir+1]*deltaT) 
                  -8.0*Epsilon[ir-1]*(Phi[ir-1]+0.5*k1[ir-1]*deltaT)     +Epsilon[ir-2]*(Phi[ir-2]+0.5*k1[ir-2]*deltaT))/(12.0*deltaR*r2[ir]);
        }
        j2[nR-2] =     Gamma[nR-2]*( Pi[nR-2]+0.5*l1[nR-2]*deltaT);
        k2[nR-2] = 0;
        //k2[nR-2] =   (-Gamma[nR-5]*( Pi[nR-5]+0.5*l1[nR-5]*deltaT)    +6.0*Gamma[nR-4]*( Pi[nR-4]+0.5*l1[nR-4]*deltaT)
        //         -18.0*Gamma[nR-3]*( Pi[nR-3]+0.5*l1[nR-3]*deltaT)   +10.0*Gamma[nR-2]*( Pi[nR-2]+0.5*l1[nR-2]*deltaT)
        //          +3.0*Gamma[nR-1]*( Pi[nR-1]+0.5*l1[nR-1]*deltaT))/(12.0*deltaR);
        l2[nR-2] = 0;
        //l2[nR-2] = (-Epsilon[nR-5]*(Phi[nR-5]+0.5*k1[nR-5]*deltaT)  +6.0*Epsilon[nR-4]*(Phi[nR-4]+0.5*k1[nR-4]*deltaT)
        //       -18.0*Epsilon[nR-3]*(Phi[nR-3]+0.5*k1[nR-3]*deltaT) +10.0*Epsilon[nR-2]*(Phi[nR-2]+0.5*k1[nR-2]*deltaT)
        //        +3.0*Epsilon[nR-1]*(Phi[nR-1]+0.5*k1[nR-1]*deltaT))/(12.0*deltaR);
        j2[nR-1] = Gamma[nR-1]*(Pi[nR-1]+0.5*l1[nR-1]*deltaT);
        k2[nR-1] = 0;
        //k2[nR-1] = (phi[nR-1]+0.5*j1[nR-1]*deltaT)/r2[nR-1] -(Phi[nR-1]+0.5*k1[nR-1]*deltaT)/r[nR-1]
        //           -0.5*(3*(Phi[nR-1]+0.5*k1[nR-1]*deltaT) -4*(Phi[nR-2]+0.5*k1[nR-2]*deltaT) +Phi[nR-3]+0.5*k1[nR-3]*deltaT); 
        l2[nR-1] = 0;
        //l2[nR-1] = -2*Gamma[nR-1]*(phi[nR-1]+0.5*j1[nR-1]*deltaT+Gamma[nR-1]*(Pi[nR-1]+0.5*l1[nR-1]*deltaT) +
        //            0.5*r[nR-1]*(Phi[nR-1]+0.5*k1[nR-1]*deltaT))/r[nR-1] +
        //           -0.5*(3*Gamma[nR-1]*Gamma[nR-1]*(Pi[nR-1]+0.5*l1[nR-1]*deltaT) -4*Gamma[nR-2]*Gamma[nR-2]*(Pi[nR-2]+0.5*l1[nR-1]*deltaT) +
        //            Gamma[nR-3]*Gamma[nR-3]*(Pi[nR-3]+0.5*l1[nR-1]*deltaT)) +
        //           -0.5*(phi[nR-1]+0.5*j1[nR-1]*deltaT)*(3*Gamma[nR-1] -4*Gamma[nR-2] +Gamma[nR-3]) ;


        //calculate j3, k3 and l3
        j3[0] =          Gamma[0]*( Pi[0]+0.5*l2[0]*deltaT);
        k3[0] =   (-25.0*Gamma[0]*( Pi[0]+0.5*l2[0]*deltaT)   +48.0*Gamma[1]*( Pi[1]+0.5*l2[1]*deltaT)   -36.0*Gamma[2]*( Pi[2]+0.5*l2[2]*deltaT) 
                   +16.0*Gamma[3]*( Pi[3]+0.5*l2[3]*deltaT)    -3.0*Gamma[4]*( Pi[4]+0.5*l2[4]*deltaT))/(12.0*deltaR);
        l3[0] = (-25.0*Epsilon[0]*(Phi[0]+0.5*k2[0]*deltaT) +48.0*Epsilon[1]*(Phi[1]+0.5*k2[1]*deltaT) -36.0*Epsilon[2]*(Phi[2]+0.5*k2[2]*deltaT)
                 +16.0*Epsilon[3]*(Phi[3]+0.5*k2[3]*deltaT)  -3.0*Epsilon[4]*(Phi[4]+0.5*k2[4]*deltaT))/(12.0*deltaR*r2[0]);
        j3[1] =    Gamma[1]*( Pi[1]+0.5*l2[1]*deltaT);
        k3[1] =   (Gamma[4]*( Pi[4]+0.5*l2[4]*deltaT)   -6.0*Gamma[3]*( Pi[3]+0.5*l2[3]*deltaT)   +18.0*Gamma[2]*( Pi[2]+0.5*l2[2]*deltaT)
             -10.0*Gamma[1]*( Pi[1]+0.5*l2[1]*deltaT)   -3.0*Gamma[0]*( Pi[0]+0.5*l2[0]*deltaT))/(12.0*deltaR);
        l3[1] = (Epsilon[4]*(Phi[4]+0.5*k2[4]*deltaT) -6.0*Epsilon[3]*(Phi[3]+0.5*k2[3]*deltaT) +18.0*Epsilon[2]*(Phi[2]+0.5*k2[2]*deltaT)
           -10.0*Epsilon[1]*(Phi[1]+0.5*k2[1]*deltaT) -3.0*Epsilon[0]*(Phi[0]+0.5*k2[0]*deltaT))/(12.0*deltaR);
        #pragma omp parallel for
        for(int ir=2;ir<nR-2;ir++){
            j3[ir] =     Gamma[ir  ]*( Pi[ir  ]+0.5*l2[ir  ]*deltaT);
            k3[ir] =   (-Gamma[ir+2]*( Pi[ir+2]+0.5*l2[ir+2]*deltaT)   +8.0*Gamma[ir+1]*( Pi[ir+1]+0.5*l2[ir+1]*deltaT) 
                    -8.0*Gamma[ir-1]*( Pi[ir-1]+0.5*l2[ir-1]*deltaT)       +Gamma[ir-2]*( Pi[ir-2]+0.5*l2[ir-2]*deltaT))/(12.0*deltaR);
            l3[ir] = (-Epsilon[ir+2]*(Phi[ir+2]+0.5*k2[ir+2]*deltaT) +8.0*Epsilon[ir+1]*(Phi[ir+1]+0.5*k2[ir+1]*deltaT) 
                  -8.0*Epsilon[ir-1]*(Phi[ir-1]+0.5*k2[ir-1]*deltaT)     +Epsilon[ir-2]*(Phi[ir-2]+0.5*k2[ir-2]*deltaT))/(12.0*deltaR*r2[ir]);
        }
        j3[nR-2] =     Gamma[nR-2]*( Pi[nR-2]+0.5*l2[nR-2]*deltaT);
        k3[nR-2] = 0;
        //k3[nR-2] =   (-Gamma[nR-5]*( Pi[nR-5]+0.5*l2[nR-5]*deltaT)    +6.0*Gamma[nR-4]*( Pi[nR-4]+0.5*l2[nR-4]*deltaT)
        //         -18.0*Gamma[nR-3]*( Pi[nR-3]+0.5*l2[nR-3]*deltaT)   +10.0*Gamma[nR-2]*( Pi[nR-2]+0.5*l2[nR-2]*deltaT)
        //          +3.0*Gamma[nR-1]*( Pi[nR-1]+0.5*l2[nR-1]*deltaT))/(12.0*deltaR);
        l3[nR-2] = 0;
        //l3[nR-2] = (-Epsilon[nR-5]*(Phi[nR-5]+0.5*k2[nR-5]*deltaT)  +6.0*Epsilon[nR-4]*(Phi[nR-4]+0.5*k2[nR-4]*deltaT)
        //       -18.0*Epsilon[nR-3]*(Phi[nR-3]+0.5*k2[nR-3]*deltaT) +10.0*Epsilon[nR-2]*(Phi[nR-2]+0.5*k2[nR-2]*deltaT)
        //        +3.0*Epsilon[nR-1]*(Phi[nR-1]+0.5*k2[nR-1]*deltaT))/(12.0*deltaR);
        j3[nR-1] = Gamma[nR-1]*(Pi[nR-1]+0.5*l2[nR-1]*deltaT);
        k3[nR-1] = 0;
        //k3[nR-1] = (phi[nR-1]+0.5*j2[nR-1]*deltaT)/r2[nR-1] -(Phi[nR-1]+0.5*k2[nR-1]*deltaT)/r[nR-1]
        //           -0.5*(3*(Phi[nR-1]+0.5*k2[nR-1]*deltaT) -4*(Phi[nR-2]+0.5*k2[nR-2]*deltaT) +Phi[nR-3]+0.5*k2[nR-3]*deltaT); 
        l3[nR-1] = 0;
        //l3[nR-1] = -2*Gamma[nR-1]*(phi[nR-1]+0.5*j2[nR-1]*deltaT+Gamma[nR-1]*(Pi[nR-1]+0.5*l2[nR-1]*deltaT) +
        //            0.5*r[nR-1]*(Phi[nR-1]+0.5*k2[nR-1]*deltaT))/r[nR-1] +
        //           -0.5*(3*Gamma[nR-1]*Gamma[nR-1]*(Pi[nR-1]+0.5*l2[nR-1]*deltaT) -4*Gamma[nR-2]*Gamma[nR-2]*(Pi[nR-2]+0.5*l2[nR-1]*deltaT) +
        //            Gamma[nR-3]*Gamma[nR-3]*(Pi[nR-3]+0.5*l2[nR-1]*deltaT)) +
        //           -0.5*(phi[nR-1]+0.5*j2[nR-1]*deltaT)*(3*Gamma[nR-1] -4*Gamma[nR-2] +Gamma[nR-3]) ;


        //calculate j4, k4 and l4
        j4[0] = Gamma[0]*(Pi[0]+l3[0]*deltaT);
        k4[0] =   (-25.0*Gamma[0]*( Pi[0]+l3[0]*deltaT)   +48.0*Gamma[1]*( Pi[1]+l3[1]*deltaT)   -36.0*Gamma[2]*( Pi[2]+l3[2]*deltaT) 
                   +16.0*Gamma[3]*( Pi[3]+l3[3]*deltaT)    -3.0*Gamma[4]*( Pi[4]+l3[4]*deltaT))/(12.0*deltaR);
        l4[0] = (-25.0*Epsilon[0]*(Phi[0]+k3[0]*deltaT) +48.0*Epsilon[1]*(Phi[1]+k3[1]*deltaT) -36.0*Epsilon[2]*(Phi[2]+k3[2]*deltaT)
                 +16.0*Epsilon[3]*(Phi[3]+k3[3]*deltaT)  -3.0*Epsilon[4]*(Phi[4]+k3[4]*deltaT))/(12.0*deltaR*r2[0]);
        j4[1] =    Gamma[1]*( Pi[1]+l3[1]*deltaT);
        k4[1] =   (Gamma[4]*( Pi[4]+l3[4]*deltaT)   -6.0*Gamma[3]*( Pi[3]+l3[3]*deltaT)   +18.0*Gamma[2]*( Pi[2]+l3[2]*deltaT)
             -10.0*Gamma[1]*( Pi[1]+l3[1]*deltaT)   -3.0*Gamma[0]*( Pi[0]+l3[0]*deltaT))/(12.0*deltaR);
        l4[1] = (Epsilon[4]*(Phi[4]+k3[4]*deltaT) -6.0*Epsilon[3]*(Phi[3]+k3[3]*deltaT) +18.0*Epsilon[2]*(Phi[2]+k3[2]*deltaT)
           -10.0*Epsilon[1]*(Phi[1]+k3[1]*deltaT) -3.0*Epsilon[0]*(Phi[0]+k3[0]*deltaT))/(12.0*deltaR);
        #pragma omp parallel for
        for(int ir=2;ir<nR-2;ir++){
            j4[ir] =     Gamma[ir  ]*( Pi[ir  ]+l3[ir  ]*deltaT);
            k4[ir] =   (-Gamma[ir+2]*( Pi[ir+2]+l3[ir+2]*deltaT)   +8.0*Gamma[ir+1]*( Pi[ir+1]+l3[ir+1]*deltaT) 
                    -8.0*Gamma[ir-1]*( Pi[ir-1]+l3[ir-1]*deltaT)       +Gamma[ir-2]*( Pi[ir-2]+l3[ir-2]*deltaT))/(12.0*deltaR);
            l4[ir] = (-Epsilon[ir+2]*(Phi[ir+2]+k3[ir+2]*deltaT) +8.0*Epsilon[ir+1]*(Phi[ir+1]+k3[ir+1]*deltaT) 
                  -8.0*Epsilon[ir-1]*(Phi[ir-1]+k3[ir-1]*deltaT)     +Epsilon[ir-2]*(Phi[ir-2]+k3[ir-2]*deltaT))/(12.0*deltaR*r2[ir]);
        }
        j4[nR-2] =     Gamma[nR-2]*( Pi[nR-2]+l3[nR-2]*deltaT);
        k4[nR-2] = 0;
        //k4[nR-2] =   (-Gamma[nR-5]*( Pi[nR-5]+l3[nR-5]*deltaT)    +6.0*Gamma[nR-4]*( Pi[nR-4]+l3[nR-4]*deltaT)
        //         -18.0*Gamma[nR-3]*( Pi[nR-3]+l3[nR-3]*deltaT)   +10.0*Gamma[nR-2]*( Pi[nR-2]+l3[nR-2]*deltaT)
        //          +3.0*Gamma[nR-1]*( Pi[nR-1]+l3[nR-1]*deltaT))/(12.0*deltaR);
        l4[nR-2] = 0;
        //l4[nR-2] = (-Epsilon[nR-5]*(Phi[nR-5]+k3[nR-5]*deltaT)  +6.0*Epsilon[nR-4]*(Phi[nR-4]+k3[nR-4]*deltaT)
        //       -18.0*Epsilon[nR-3]*(Phi[nR-3]+k3[nR-3]*deltaT) +10.0*Epsilon[nR-2]*(Phi[nR-2]+k3[nR-2]*deltaT)
        //        +3.0*Epsilon[nR-1]*(Phi[nR-1]+k3[nR-1]*deltaT))/(12.0*deltaR);
        j4[nR-1] = Gamma[nR-1]*(Pi[nR-1]+l3[nR-1]*deltaT);
        k4[nR-1] = 0;
        //k4[nR-1] = (phi[nR-1]+j3[nR-1]*deltaT)/r2[nR-1] -(Phi[nR-1]+k3[nR-1]*deltaT)/r[nR-1]
        //           -0.5*(3*(Phi[nR-1]+k3[nR-1]*deltaT) -4*(Phi[nR-2]+k3[nR-2]*deltaT) +Phi[nR-3]+k3[nR-3]*deltaT); 
        l4[nR-1] = 0;
        //l4[nR-1] = -2*Gamma[nR-1]*(phi[nR-1]+j3[nR-1]*deltaT+Gamma[nR-1]*(Pi[nR-1]+l3[nR-1]*deltaT) +
        //            0.5*r[nR-1]*(Phi[nR-1]+k3[nR-1]*deltaT))/r[nR-1] +
        //           -0.5*(3*Gamma[nR-1]*Gamma[nR-1]*(Pi[nR-1]+l3[nR-1]*deltaT) -4*Gamma[nR-2]*Gamma[nR-2]*(Pi[nR-2]+l3[nR-1]*deltaT) +
        //            Gamma[nR-3]*Gamma[nR-3]*(Pi[nR-3]+l3[nR-1]*deltaT)) +
        //           -0.5*(phi[nR-1]+j3[nR-1]*deltaT)*(3*Gamma[nR-1] -4*Gamma[nR-2] +Gamma[nR-3]) ;


        //Calculate phi, Phi and Pi on next step
        #pragma omp parallel for
        for(int ir=0;ir<nR;ir++){
            phi[ir] += (j1[ir]+2.0*(j2[ir]+j3[ir])+j4[ir])*deltaT/6.0;
            Phi[ir] += (k1[ir]+2.0*(k2[ir]+k3[ir])+k4[ir])*deltaT/6.0;
            Pi[ir]  += (l1[ir]+2.0*(l2[ir]+l3[ir])+l4[ir])*deltaT/6.0;
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

void print_mass(double *MassHist,struct tm time_data,int print_iterations){
    char fileName[50];
    snprintf(fileName, sizeof(fileName), "Mass_%02d%02d%02d.dat", time_data.tm_hour, time_data.tm_min, time_data.tm_sec);
    FILE* data = fopen(fileName,"w");
    //Print M
    for(int i=0;i<print_iterations-1;i++){
        fprintf(data,"%lf,",MassHist[i]);
    }
    fprintf(data,"%lf",MassHist[print_iterations-1]);
    fclose(data);
}

void print_data(double** hist,int fType,double* model_parameters,int iterations,int maxR,double deltaR,int nP,time_t totalTime){
    int print_iterations, printR;
    if(SAVE_MODE == 0){
        print_iterations = iterations/SAVE_ITERATION;
        printR = (maxR/deltaR)/SAVE_RES;
    }
    if(SAVE_MODE == 1){
        print_iterations = iterations;
        printR = maxR/deltaR;
    }
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
    double q  = model_parameters[3];
    //Add time to filename
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    if(MASS) print_mass(hist[6],tm,print_iterations);
    char fileName[50];
    snprintf(fileName, sizeof(fileName), "Output_%02d%02d%02d.dat", tm.tm_hour, tm.tm_min, tm.tm_sec);
    FILE* data = fopen(fileName,"w");
    //Print all parameters 
         if(METRIC == 0) fprintf(data,"Metric: Minkowski\n");
    else if(METRIC == 1) fprintf(data,"Metric: Choptuik\n");
    else if(METRIC == 2) fprintf(data,"Metric: Modified Choptuik\n" );
    else if(METRIC == 3) fprintf(data,"Metric: Choptuik Euler\n");
    else if(METRIC == 4) fprintf(data,"Metric: Modified Choptuik Euler\n" );
         if(fType == 0) fprintf(data,"Function type: Hyperbolic Tan\n");
    else if(fType == 1) fprintf(data,"Function type: Exponential\n");
    else if(fType == 2) fprintf(data,"Function type: Constant\n");
    else if(fType == 3) fprintf(data,"Function type: Linear\n");
    fprintf(data,"p0: %lf\n",p0);
    fprintf(data,"r0: %lf\n",r0);
    fprintf(data,"d: %lf\n",d);
    fprintf(data,"q: %lf\n",q);
    fprintf(data,"R step size: %lf\n",deltaR);
    fprintf(data,"Maximum R: %d\n",maxR);
    fprintf(data,"Iterations: %d\n",iterations);
    fprintf(data,"Number of threads: %d\n",nP);
    fprintf(data,"Total simulation time: %ld\n",totalTime);
    //Print R
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(data,"%lf,",hist[0][ir]);
    }
    fprintf(data,"%lf\n",hist[0][printR-1]);
    //Print phi
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[1][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[1][i*printR+printR-1]);
    }
    //Print Phi
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[2][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[2][i*printR+printR-1]);
    }
    //Print Pi
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[3][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[3][i*printR+printR-1]);
    }
    //Print a
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[4][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[4][i*printR+printR-1]);
    }
    //Print alpha
    for(int i=0;i<print_iterations-1;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[5][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[5][i*printR+printR-1]);
    }
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(data,"%lf,",hist[5][(print_iterations-1)*printR+ir]);
    }
    fprintf(data,"%lf",hist[5][print_iterations*printR-1]);
    fclose(data);
}


int main(int argc, char* argv[]){
    double* model_parameters = malloc(sizeof(double)*4);
    double** initial_conditions;
    double** hist;
    double* r;
    double* phi;
    double* Phi;
    double* Pi;
    //Define simulation parameters
    int fType = 0;
    if((argc>1) && atoi(argv[1])) fType = atoi(argv[1]);
    double p0 = 0.001;
    if((argc>2) && atof(argv[2])) p0 = atof(argv[2]);
    double r0 = 20.;
    if((argc>3) && atof(argv[3])) r0 = atof(argv[3]);
    double d = 3.;
    if((argc>4) && atof(argv[4])) d  = atof(argv[4]);
    double q = 2.;
    if((argc>5) && atof(argv[5])) q  = atof(argv[5]);
    model_parameters[0] = p0;
    model_parameters[1] = r0;
    model_parameters[2] = d;
    model_parameters[3] = q;
    
    //Define simulation limits
    double deltaR = 0.001;
    if((argc>6) && atof(argv[6])) deltaR = atof(argv[6]);
    int maxR = 70;
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
    hist = iteration(r,phi,Phi,Pi,deltaR,maxR,iterations,SAVE_ITERATION);
    time_t finalTime = time(NULL);
    int nP = omp_get_max_threads();
    time_t timeDelta = (finalTime-initTime);

    //Print simulation history to a file
    if(SAVE_MODE == 0)
        print_data(hist,fType,model_parameters,iterations,maxR,deltaR,nP,timeDelta);
    else if(SAVE_MODE == 1){
        print_data(hist,fType,model_parameters,iterations-FIRST_ITERATION,maxR-MIN_R,deltaR,nP,timeDelta);
        }
    printf("Finished, total time: %lds\n", timeDelta);
}




