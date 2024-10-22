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

#define MASS 0
#define METRIC 1 // 0 = minkowski; 1 = choptuik; 2 = modified choptuik
#define SAVE_MODE 1 // 0 = save uniformly on every SAVE_RES and SAVE_ITERATION ; 1 = save all points after FIRST_ITERATION and with r > MIN_R
#define SAVE_RES 50
#define SAVE_ITERATION 10
#define FIRST_ITERATION 7700 //77950
#define MIN_R 0
#define ITERATIONS 78500
#define EPSILON 0.0 // Default Kreiss-Oliger dampening
#define PI 3.141592653
#define E  2.718281828
#define _6 0.166666666 // 1/6

double leftmost_D1(double *f,int idx,double deltaR){
    return (-25.0*f[idx] +48.0*f[idx+1] -36.0*f[idx+2] +16.0*f[idx+3] -3.0*f[idx+4])/(12.0*deltaR);
}
double leftmid_D1(double *f,int idx,double deltaR){
    return (-3.0*f[idx-1] -10.0*f[idx] +18.0*f[idx+1] -6.0*f[idx+2] +1.0*f[idx+3])/(12.0*deltaR);
}
double centered_D1(double *f,int idx,double deltaR){
    return (f[idx-2] -8.0*f[idx-1] +8.0*f[idx+1] -f[idx+2])/(12.0*deltaR);
}
double rightmid_D1(double *f,int idx,double deltaR){
    return (-1.0*f[idx-3] +6.0*f[idx-2] -18.0*f[idx-1] +10.0*f[idx] +3.0*f[idx+1])/(12.0*deltaR);
}
double rightmost_D1(double *f,int idx,double deltaR){
    return (3.0*f[idx-4] -16.0*f[idx-3] +36.0*f[idx-2] -48.0*f[idx-1] +25.0*f[idx])/(12.0*deltaR);
}
double centered_r0_odd_D1(double *f,double deltaR){
    return (16.0*f[1] -2.0*f[2])/(12.0*deltaR);
}
double centered_r0_even_D1(double *f,double deltaR){
    return 0.0;
}
double centered_r1_odd_D1(double *f,double deltaR){
    return (-f[1] -8.0*f[0] +8.0*f[2] -f[3])/(12.0*deltaR);
}
double centered_r1_even_D1(double *f,double deltaR){
    return ( f[1] -8.0*f[0] +8.0*f[2] -f[3])/(12.0*deltaR);
}

double** initialize_field(int fType,double* model_parameters,double deltaR,int maxR){
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
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
            phi[i] = p0*pow(r[i],3)*pow(E,-pow((r[i]-r0)/d,2));
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

void metric_iteration1(double* Beta, double* Beta1_2, double* a, double* alpha, double* r, int nR, double deltaR){
    
    int m1, n1, m2, n2, m3, n3, m4, n4;
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
}

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
    double *Zeta = malloc(sizeof(double)*nR);
    double *Phi_ko = malloc(sizeof(double)*nR);
    double *Pi_ko = malloc(sizeof(double)*nR);
    double ko_c = pow(-1,3)*epsilon/5.;
    double *r2 = malloc(sizeof(double)*nR);

    double *Chi = malloc(sizeof(double)*5);
    double *j_n = malloc(sizeof(double)*nR);
    double *k_n = malloc(sizeof(double)*nR);
    double *l_n = malloc(sizeof(double)*nR);
    double *j_sum = malloc(sizeof(double)*nR);
    double *k_sum = malloc(sizeof(double)*nR);
    double *l_sum = malloc(sizeof(double)*nR);
    double  rk[4] = {1.0,2.0,2.0,1.0};
    double _rk[4] = {1.0,0.5,0.5,1.0};

    double m1, n1, m2, n2, m3, n3, m4, n4;
    double *m, *n;
    if(METRIC==5){
        m = malloc(sizeof(double)*4*nR);
        n = malloc(sizeof(double)*4);
    }
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
    if(METRIC == 0){
        for(int i=1;i<nR;i++){
            //a[i] = 1.0 + 0.0001*pow(r[i],3)*pow(E,-pow((r[i]-10.0)/10,2));
            //alpha[i] = pow(E,r[i]/100);
        }
    }
    for(int ir=0;ir<nR;ir++){
        Gamma[ir] = alpha[ir]/a[ir];
        r2[ir] = r[ir]*r[ir];
        Epsilon[ir] = r2[ir]*Gamma[ir];
        Zeta[ir] = r[ir]/(a[ir]*a[ir]);
    }

    for(int i=0;i<iterations;i++){
        //printf("Iteration %d\n",i);
        //If on minkowski metric, do not solve a and alpha

        //Useful auxiliar variable
        //#pragma omp parallel for
        for(int ir=0;ir<nR;ir++)
            Beta[ir] = 2.0*PI*r[ir]*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
        //If on choptuik metric, solve a and alpha with known values of PHI and PI using RK4
        if(METRIC == 1){
            //#pragma omp parallel for
            for(int ir=0;ir<nR-1;ir++){
                Beta1_2[ir] = 0.5*PI*(r[ir]+0.5*deltaR)*((Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])+(Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1]));
                //Beta1_2[ir] = 2.0*PI*(r[ir]+deltaR/2)*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
            }
            
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
            //Define auxiliar variables
            //#pragma omp parallel for
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
        //If on choptuik metric, solve a and alpha with known values of PHI and PI using RK4
        if(METRIC == 4){
            //#pragma omp parallel for
            for(int ir=0;ir<nR-1;ir++){
                Beta[ir] = 2.0*PI*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
                //Beta1_2[ir] = 2.0*PI*(r[ir]+deltaR/2)*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
            }
            
            for(int ir=0;ir<nR-1;ir++){
                //calculate m1 and n1
                m1 = deltaR*    a[ir]*(r[ir]*Beta[ir]-0.5*(a[ir]*a[ir]-1)/r[ir]);
                n1 = deltaR*alpha[ir]*(r[ir]*Beta[ir]+0.5*(a[ir]*a[ir]-1)/r[ir]);
                //calculate m2 and n2
                m2 = deltaR*    (a[ir]+0.5*m1)*((r[ir]+deltaR)*Beta[ir]-0.5*((a[ir]+0.5*m1)*(a[ir]+0.5*m1)-1)/(r[ir]+0.5*deltaR));
                n2 = deltaR*(alpha[ir]+0.5*n1)*((r[ir]+deltaR)*Beta[ir]+0.5*((a[ir]+0.5*m1)*(a[ir]+0.5*m1)-1)/(r[ir]+0.5*deltaR));
                //calculate m3 and n3
                m3 = deltaR*    (a[ir]+0.5*m2)*((r[ir]+deltaR)*Beta[ir]-0.5*((a[ir]+0.5*m2)*(a[ir]+0.5*m2)-1)/(r[ir]+0.5*deltaR));
                n3 = deltaR*(alpha[ir]+0.5*n2)*((r[ir]+deltaR)*Beta[ir]+0.5*((a[ir]+0.5*m2)*(a[ir]+0.5*m2)-1)/(r[ir]+0.5*deltaR));
                //calculate m4 and n4
                m4 = deltaR*    (a[ir]+m3)*(r[ir+1]*Beta[ir]-0.5*((a[ir]+m3)*(a[ir]+m3)-1)/r[ir+1]);
                n4 = deltaR*(alpha[ir]+n3)*(r[ir+1]*Beta[ir]+0.5*((a[ir]+m3)*(a[ir]+m3)-1)/r[ir+1]);
                //Calculate next step for a and alpha
                a[ir+1]     = a[ir]     +(m1 +2.0*(m2+m3) +m4)/6.0;
                alpha[ir+1] = alpha[ir] +(n1 +2.0*(n2+n3) +n4)/6.0;
            }
            //Define auxiliar variables
            //#pragma omp parallel for
            for(int ir=0;ir<nR;ir++){
                Gamma[ir] = alpha[ir]/a[ir];
                Epsilon[ir] = r2[ir]*Gamma[ir];
                Zeta[ir] = r[ir]/(a[ir]*a[ir]);
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
                for(int ir=0;ir<saveNR;ir++){
                    Xhistory[save_iter*saveNR + ir] =   Phi[noSaveNR + ir];
                    Yhistory[save_iter*saveNR + ir] =    Pi[noSaveNR + ir];
                    Fhistory[save_iter*saveNR + ir] =   phi[noSaveNR + ir];
                    Ahistory[save_iter*saveNR + ir] =     a[noSaveNR + ir];
                    Bhistory[save_iter*saveNR + ir] = alpha[noSaveNR + ir];
                }
            }
        }

        //Calculate Kreiss-Oliger dissipation of 6th order for next step
        Phi_ko[0] = (14.71804060*Phi[0]-99.3467741*Phi[1]+287.0017918*Phi[2]-459.938769*Phi[3]+441.5412182*Phi[4]-253.8862004*Phi[5]+80.9492233*Phi[6]-11.03853045*Phi[7])/3.67951015;
        Pi_ko[ 0] = (14.71804060* Pi[0]-99.3467741* Pi[1]+287.0017918* Pi[2]-459.938769* Pi[3]+441.5412182* Pi[4]-253.8862004* Pi[5]+80.9492233* Pi[6]-11.03853045* Pi[7])/3.67951015;
        Phi_ko[1] = 3*Phi[0]-20*Phi[1]+57*Phi[2]-90*Phi[3]+85*Phi[4]-48*Phi[5]+15*Phi[6]-2*Phi[7];
        Pi_ko[ 1] = 3* Pi[0]-20* Pi[1]+57* Pi[2]-90* Pi[3]+85* Pi[4]-48* Pi[5]+15* Pi[6]-2* Pi[7];
        Phi_ko[2] = 2*Phi[0]-13*Phi[1]+36*Phi[2]-55*Phi[3]+50*Phi[4]-27*Phi[5]+8*Phi[6]-1*Phi[7];
        Pi_ko[ 2] = 2* Pi[0]-13* Pi[1]+36* Pi[2]-55* Pi[3]+50* Pi[4]-27* Pi[5]+8* Pi[6]-1* Pi[7];
        //Phi_ko[0] = Phi[0]-6.0*Phi[1]+15.0*Phi[2]-20.0*Phi[3]+15.0*Phi[4]-6.0*Phi[5]+Phi[6];
        //Pi_ko[ 0] =  Pi[0]-6.0* Pi[1]+15.0* Pi[2]-20.0* Pi[3]+15.0* Pi[4]-6.0* Pi[5]+ Pi[6];
        //Phi_ko[1] = Phi[0]-6.0*Phi[1]+15.0*Phi[2]-20.0*Phi[3]+15.0*Phi[4]-6.0*Phi[5]+Phi[6];
        //Pi_ko[ 1] =  Pi[0]-6.0* Pi[1]+15.0* Pi[2]-20.0* Pi[3]+15.0* Pi[4]-6.0* Pi[5]+ Pi[6];
        //Phi_ko[2] = Phi[0]-6.0*Phi[1]+15.0*Phi[2]-20.0*Phi[3]+15.0*Phi[4]-6.0*Phi[5]+Phi[6];
        //Pi_ko[ 2] =  Pi[0]-6.0* Pi[1]+15.0* Pi[2]-20.0* Pi[3]+15.0* Pi[4]-6.0* Pi[5]+ Pi[6];
        //Phi_ko[0] = -20.0*Phi[0]+30.0*Phi[1]-12.0*Phi[2] +2.0*Phi[3];
        //Pi_ko[ 0] = -20.0* Pi[0]+30.0* Pi[1]-12.0* Pi[2] +2.0* Pi[3];
        //Phi_ko[1] =  15.0*Phi[0]-26.0*Phi[1]+16.0*Phi[2] -6.0*Phi[3]    +Phi[4];
        //Pi_ko[ 1] =  15.0* Pi[0]-26.0* Pi[1]+16.0* Pi[2] -6.0* Pi[3]    + Pi[4];
        //Phi_ko[2] =  -6.0*Phi[0]+16.0*Phi[1]-20.0*Phi[2]+15.0*Phi[3]-6.0*Phi[4]+Phi[5];
        //Pi_ko[ 2] =  -6.0* Pi[0]+16.0* Pi[1]-20.0* Pi[2]+15.0* Pi[3]-6.0* Pi[4]+ Pi[5];
        for(int ir=3;ir<nR-3;ir++){
            Phi_ko[ir] = Phi[ir-3]-6.0*Phi[ir-2]+15.0*Phi[ir-1]-20.0*Phi[ir]+15.0*Phi[ir+1]-6.0*Phi[ir+2]+Phi[ir+3];
            Pi_ko[ir]  =  Pi[ir-3]-6.0* Pi[ir-2]+15.0* Pi[ir-1]-20.0* Pi[ir]+15.0* Pi[ir+1]-6.0* Pi[ir+2]+ Pi[ir+3];
        }
        Phi_ko[nR-3] = Phi[nR-6]-6.0*Phi[nR-5]+15.0*Phi[nR-4]-20.0*Phi[nR-3]+15.0*Phi[nR-2]-6.0*Phi[nR-1];
        Pi_ko[ nR-3] =  Pi[nR-6]-6.0* Pi[nR-5]+15.0* Pi[nR-4]-20.0* Pi[nR-3]+15.0* Pi[nR-2]-6.0* Pi[nR-1];
        Phi_ko[nR-2] = Phi[nR-5]-6.0*Phi[nR-4]+15.0*Phi[nR-3]-20.0*Phi[nR-2]+15.0*Phi[nR-1];
        Pi_ko[ nR-2] =  Pi[nR-5]-6.0* Pi[nR-4]+15.0* Pi[nR-3]-20.0* Pi[nR-2]+15.0* Pi[nR-1];
        Phi_ko[nR-1] = Phi[nR-4]-6.0*Phi[nR-3]+15.0*Phi[nR-2]-20.0*Phi[nR-1];
        Pi_ko[ nR-1] =  Pi[nR-4]-6.0* Pi[nR-3]+15.0* Pi[nR-2]-20.0* Pi[nR-1];

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
                Gamma[ir]   =        alpha[ir]*( Pi[ir]+_rk[n]*deltaT*l_n[ir])/a[ir];
                Epsilon[ir] = r2[ir]*alpha[ir]*(Phi[ir]+_rk[n]*deltaT*k_n[ir])/a[ir];
            }
            for(int ir=0;ir<5;ir++)
                Chi[ir] = alpha[ir]*(Phi[ir]+_rk[n]*deltaT*k_n[ir])/a[ir];

            j_n[0] = Gamma[0];
            //k_n[0] = leftmost_D1(Gamma,0,deltaR);
            //l_n[0] = leftmost_D1(Epsilon,0,deltaR)/r2[0];
            k_n[0] = 0;
            l_n[0] = centered_r0_odd_D1(Chi,deltaR); //Chi is odd because Phi is odd
            j_n[1] = Gamma[1];
            //k_n[1] = leftmid_D1(Gamma,1,deltaR);
            //l_n[1] = leftmid_D1(Epsilon,1,deltaR)/r2[1];
            k_n[1] = centered_r1_even_D1(Gamma,deltaR); //Gamma is even because Pi es even
            l_n[1] = centered_r1_odd_D1(Epsilon,deltaR)/r2[1]; //Epsilon is odd because Phi is odd
            //#pragma omp parallel for
            for(int ir=2;ir<nR-2;ir++){
                j_n[ir] = Gamma[ir];
                k_n[ir] = centered_D1(Gamma,ir,deltaR);
                l_n[ir] = centered_D1(Epsilon,ir,deltaR)/r2[ir];
            }
            j_n[nR-2] = Gamma[nR-2];
            k_n[nR-2] = 0;
            l_n[nR-2] = 0;
            j_n[nR-1] = Gamma[nR-1];
            k_n[nR-1] = 0;
            l_n[nR-1] = 0;

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
            Phi[ir] += _6*deltaT*k_sum[ir] -ko_c*Phi_ko[ir];
            Pi[ir]  += _6*deltaT*l_sum[ir] -ko_c* Pi_ko[ir];
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

double** iteration2(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration,double epsilon){

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
    double ko_c = pow(-1,3)*epsilon/5.;
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
                        mass += Beta[ir]*r[ir]/(a[ir]*a[ir])*deltaR;
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
                        mass += Beta[ir]*r[ir]/(a[ir]*a[ir])*deltaR;
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
            for(int ir=0;ir<nR;ir++){
                Gamma[ir]   =        alpha[ir]*( Pi_rk[ir])/a[ir];
                Epsilon[ir] = r2[ir]*alpha[ir]*(Phi_rk[ir])/a[ir];
            }
            for(int ir=0;ir<5;ir++){
                rPhi[ir] = r[nR-5 +ir]*Phi_rk[nR-5 +ir];
                 rPi[ir] = r[nR-5 +ir]* Pi_rk[nR-5 +ir];
            }



            j_n[0] = Gamma[0];
            k_n[0] = 0;
            l_n[0] = leftmost_D1(Phi_rk,0,deltaR);
            //l_n[0] = centered_r0_odd_D1(Phi_rk,deltaR); //Phi is odd
            j_n[1] = Gamma[1];
            k_n[1] = leftmid_D1(Gamma,1,deltaR);
            l_n[1] = leftmid_D1(Epsilon,1,deltaR)/r2[1];
            //k_n[1] = centered_r1_even_D1(Gamma,deltaR); //Gamma is even because Pi es even
            //l_n[1] = centered_r1_odd_D1(Epsilon,deltaR)/r2[1]; //Epsilon is odd because Phi is odd
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
            Phi[ir] += _6*deltaT*k_sum[ir];//-ko_c*Phi_ko[ir];
            Pi[ir]  += _6*deltaT*l_sum[ir];//-ko_c* Pi_ko[ir];
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

void print_data(double** hist,int fType,double* model_parameters,int iterations,int maxR,double deltaR,int nP,time_t totalTime,double epsilon){
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
    else if(METRIC == 2) fprintf(data,"Metric: Modified Choptuik\n");
    else if(METRIC == 3) fprintf(data,"Metric: Choptuik Euler\n");
    else if(METRIC == 4) fprintf(data,"Metric: Choptuik Alt\n");
    else                 fprintf(data,"Metric: Other\n");
         if(fType == 0) fprintf(data,"Function type: Hyperbolic Tan\n");
    else if(fType == 1) fprintf(data,"Function type: Exponential\n");
    else if(fType == 2) fprintf(data,"Function type: Constant\n");
    else if(fType == 3) fprintf(data,"Function type: Linear\n");
    fprintf(data,"p0: %lf\n",p0);
    fprintf(data,"r0: %lf\n",r0);
    fprintf(data,"d: %lf\n",d);
    fprintf(data,"Kreiss-Oliger Coefficient: %lf\n",epsilon);
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
    hist = iteration2(r,phi,Phi,Pi,deltaR,maxR,iterations,SAVE_ITERATION,epsilon);
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




