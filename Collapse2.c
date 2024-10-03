#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define MASS 0
#define METRIC 1 // 0 = minkowski; 1 = choptuik; 2 = modified choptuik
#define SAVE_RES 500
#define SAVE_ITERATION 100
#define FIRST_ITERATION 77950
#define MIN_R 50
#define ITERATIONS 78100
#define EPSILON 0.0 // Default Kreiss-Oliger dampening, must be smaller than 0.15625
#define PI 3.141592653
#define E  2.718281828

double** initialize_field(double* model_parameters,double deltaR,int maxR){
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

    for(int i=0;i<nR;i++)
        phi[i] = p0*pow(r[i],3)*pow(E,-pow((r[i]-r0)/d,2));

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

void metric_iterationRK1(double *r,double *Phi,double *Pi,double *a, double *alpha,double deltaR,int nR){
    double m1, n1;

    for(int ir=0;ir<nR-1;ir++){
        m1 = deltaR*    a[ir]*(2.0*PI*r[ir]*(Phi[ir]*Phi[ir]+Pi[ir]*Pi[ir]) -0.5*(a[ir]*a[ir]-1)/r[ir]);
        n1 = deltaR*alpha[ir]*(2.0*PI*r[ir]*(Phi[ir]*Phi[ir]+Pi[ir]*Pi[ir]) +0.5*(a[ir]*a[ir]-1)/r[ir]);

            a[ir+1] += m1;
        alpha[ir+1] += n1;
    }
}
void metric_iterationRK4(double *r,double *Phi,double *Pi,double *a, double *alpha,double deltaR,int nR){
    double m1, n1, m2, n2, m3, n3, m4, n4;

    for(int ir=0;ir<nR-1;ir++){
        m1 = deltaR*    a[ir]*(2.0*PI*r[ir]*(Phi[ir]*Phi[ir]+Pi[ir]*Pi[ir]) -0.5*(a[ir]*a[ir]-1)/r[ir]);
        n1 = deltaR*alpha[ir]*(2.0*PI*r[ir]*(Phi[ir]*Phi[ir]+Pi[ir]*Pi[ir]) +0.5*(a[ir]*a[ir]-1)/r[ir]);

        m2 = deltaR*(    a[ir]+0.5*m1)*(0.5*PI*(r[ir]+0.5*deltaR)*((Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1])+(Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])) -0.5*((a[ir]+0.5*m1)*(a[ir]+0.5*m1)-1)/(r[ir]+0.5*deltaR));
        n2 = deltaR*(alpha[ir]+0.5*n1)*(0.5*PI*(r[ir]+0.5*deltaR)*((Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1])+(Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])) +0.5*((a[ir]+0.5*m1)*(a[ir]+0.5*m1)-1)/(r[ir]+0.5*deltaR));

        m3 = deltaR*(    a[ir]+0.5*m2)*(0.5*PI*(r[ir]+0.5*deltaR)*((Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1])+(Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])) -0.5*((a[ir]+0.5*m2)*(a[ir]+0.5*m2)-1)/(r[ir]+0.5*deltaR));
        n3 = deltaR*(alpha[ir]+0.5*n2)*(0.5*PI*(r[ir]+0.5*deltaR)*((Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1])+(Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])) +0.5*((a[ir]+0.5*m2)*(a[ir]+0.5*m2)-1)/(r[ir]+0.5*deltaR));

        m4 = deltaR*(    a[ir]+m3)*(2.0*PI*(r[ir]+deltaR)*(Phi[ir+1]*Phi[ir+1]+Pi[ir+1]*Pi[ir+1]) -0.5*((a[ir]+m3)*(a[ir]+m3)-1)/(r[ir]+deltaR));
        n4 = deltaR*(alpha[ir]+n3)*(2.0*PI*(r[ir]+deltaR)*(Phi[ir+1]*Phi[ir+1]+Pi[ir+1]*Pi[ir+1]) +0.5*((a[ir]+m3)*(a[ir]+m3)-1)/(r[ir]+deltaR));

            a[ir+1] += (m1 +2.0*m2 +2.0*m3 +m4)/6.0;
        alpha[ir+1] += (n1 +2.0*n2 +2.0*n3 +n4)/6.0;
    }
}

double **iteration(double *r,double *phi,double *Phi,double *Pi,double deltaR,int nR,int iterations){
    double deltaT = deltaR/5.0;
    double *rHist = malloc(sizeof(double)*(nR/SAVE_RES));
    double *phiHist = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/SAVE_ITERATION));
    double *PhiHist = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/SAVE_ITERATION));
    double *PiHist = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/SAVE_ITERATION));
    double *aHist = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/SAVE_ITERATION));
    double *alphaHist = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/SAVE_ITERATION));
    double *a = malloc(sizeof(double)*nR);
    double *alpha = malloc(sizeof(double)*nR);
    double *j1 = malloc(sizeof(double)*nR);
    double *k1 = malloc(sizeof(double)*nR);
    double *l1 = malloc(sizeof(double)*nR);
    double *j2 = malloc(sizeof(double)*nR);
    double *k2 = malloc(sizeof(double)*nR);
    double *l2 = malloc(sizeof(double)*nR);
    double *j3 = malloc(sizeof(double)*nR);
    double *k3 = malloc(sizeof(double)*nR);
    double *l3 = malloc(sizeof(double)*nR);
    double *j4 = malloc(sizeof(double)*nR);
    double *k4 = malloc(sizeof(double)*nR);
    double *l4 = malloc(sizeof(double)*nR);
    double **history = malloc(sizeof(double*)*6);
    int save_count = SAVE_ITERATION;



    for(int ir=0;ir<nR;ir++){
        a[ir] = 1.0;
        alpha[ir] = 1.0;
    }
    for(int ir=0;ir<(nR/SAVE_RES);ir++){
        rHist[ir] = r[ir*SAVE_RES];
    }
    for(int i;i<iterations;i++){
        //Iterate the metric
        metric_iterationRK4(r,Phi,Pi,a,alpha,deltaR,nR);

        //Save values of Phi, Pi, phi, a and alpha
        if(save_count == SAVE_ITERATION){
            //printf("iteration %d\n",i);
            //#pragma omp parallel for
            for(int ir=0;ir<(nR/SAVE_RES);ir++){
                  phiHist[(i/SAVE_ITERATION)*(nR/SAVE_RES)+(ir)] =   phi[ir*SAVE_RES];
                  PhiHist[(i/SAVE_ITERATION)*(nR/SAVE_RES)+(ir)] =   Phi[ir*SAVE_RES];
                   PiHist[(i/SAVE_ITERATION)*(nR/SAVE_RES)+(ir)] =    Pi[ir*SAVE_RES];
                    aHist[(i/SAVE_ITERATION)*(nR/SAVE_RES)+(ir)] =     a[ir*SAVE_RES];
                alphaHist[(i/SAVE_ITERATION)*(nR/SAVE_RES)+(ir)] = alpha[ir*SAVE_RES];
            }
            save_count=0;
        }
        save_count+=1;

        //Iterate phi, Phi and Pi
        j1[0] = deltaT*Pi[0]*alpha[0]/a[0];
        k1[0] = deltaT*(-25.0*alpha[0]*Pi[0]/a[0] +48.0*alpha[1]*Pi[1]/a[1] 
                        -36.0*alpha[2]*Pi[2]/a[2] +16.0*alpha[3]*Pi[3]/a[3]
                         -3.0*alpha[4]*Pi[4]/a[4])/(12.0*deltaR);
        l1[0] = deltaT*(-25.0*r[0]*r[0]*alpha[0]*Phi[0]/a[0] +48.0*r[1]*r[1]*alpha[1]*Phi[1]/a[1]
                        -36.0*r[2]*r[2]*alpha[2]*Phi[2]/a[2] +16.0*r[3]*r[3]*alpha[3]*Phi[3]/a[3]
                         -3.0*r[4]*r[4]*alpha[4]*Phi[4]/a[4])/(12.0*deltaR*r[0]*r[0]);
        j1[1] = deltaT*Pi[1]*alpha[1]/a[1];
        k1[1] = deltaT*(-3.0*alpha[0]*Pi[0]/a[0] -10.0*alpha[1]*Pi[1]/a[1]
                       +18.0*alpha[2]*Pi[2]/a[2]  -6.0*alpha[3]*Pi[3]/a[3]
                            +alpha[4]*Pi[4]/a[4])/(12.0*deltaR);
        l1[1] = deltaT*(-3.0*r[0]*r[0]*alpha[0]*Pi[0]/a[0] -10.0*r[1]*r[1]*alpha[1]*Pi[1]/a[1]
                       +18.0*r[2]*r[2]*alpha[2]*Pi[2]/a[2]  -6.0*r[3]*r[3]*alpha[3]*Pi[3]/a[3]
                            +r[4]*r[4]*alpha[4]*Pi[4]/a[4])/(12.0*deltaR*r[1]*r[1]);
        for(int ir=2;ir<nR-2;ir++){
            j1[ir] = deltaT*Pi[ir]*alpha[ir]/a[ir];
            k1[ir] = deltaT*(alpha[ir-2]*Pi[ir-2]/a[ir-2] -8.0*alpha[ir-1]*Pi[ir-1]/a[ir-1]
                        +8.0*alpha[ir+1]*Pi[ir+1]/a[ir+1]     -alpha[ir+2]*Pi[ir+2]/a[ir+2])/(12.0*deltaR);
            l1[ir] = deltaT*(r[ir-2]*r[ir-2]*alpha[ir-2]*Phi[ir-2]/a[ir-2] -8.0*r[ir-1]*r[ir-1]*alpha[ir-1]*Phi[ir-1]/a[ir-1]
                        +8.0*r[ir+1]*r[ir+1]*alpha[ir+1]*Phi[ir+1]/a[ir+1]     -r[ir+2]*r[ir+2]*alpha[ir+2]*Phi[ir+2]/a[ir+2])/(12.0*deltaR*r[ir]*r[ir]);
        }
        j1[nR-2] = 0;
        k1[nR-2] = 0;
        l1[nR-2] = 0;
        j1[nR-1] = 0;
        k1[nR-1] = 0;
        l1[nR-1] = 0;

        j2[0] = deltaT*(Pi[0]+0.5*l1[0])*alpha[0]/a[0];
        k2[0] = deltaT*(-25.0*alpha[0]*(Pi[0]+0.5*l1[0])/a[0] +48.0*alpha[1]*(Pi[1]+0.5*l1[1])/a[1] 
                        -36.0*alpha[2]*(Pi[2]+0.5*l1[2])/a[2] +16.0*alpha[3]*(Pi[3]+0.5*l1[3])/a[3]
                         -3.0*alpha[4]*(Pi[4]+0.5*l1[4])/a[4])/(12.0*deltaR);
        l2[0] = deltaT*(-25.0*r[0]*r[0]*alpha[0]*(Phi[0]+0.5*k1[0])/a[0] +48.0*r[1]*r[1]*alpha[1]*(Phi[1]+0.5*k1[1])/a[1]
                        -36.0*r[2]*r[2]*alpha[2]*(Phi[2]+0.5*k1[2])/a[2] +16.0*r[3]*r[3]*alpha[3]*(Phi[3]+0.5*k1[3])/a[3]
                         -3.0*r[4]*r[4]*alpha[4]*(Phi[4]+0.5*k1[4])/a[4])/(12.0*deltaR*r[0]*r[0]);
        j2[1] = deltaT*(Pi[1]+0.5*l1[1])*alpha[1]/a[1];
        k2[1] = deltaT*(-3.0*alpha[0]*(Pi[0]+0.5*l1[0])/a[0] -10.0*alpha[1]*(Pi[1]+0.5*l1[1])/a[1]
                       +18.0*alpha[2]*(Pi[2]+0.5*l1[2])/a[2]  -6.0*alpha[3]*(Pi[3]+0.5*l1[3])/a[3]
                            +alpha[4]*(Pi[4]+0.5*l1[4])/a[4])/(12.0*deltaR);
        l2[1] = deltaT*(-3.0*r[0]*r[0]*alpha[0]*(Phi[0]+0.5*k1[0])/a[0] -10.0*r[1]*r[1]*alpha[1]*(Phi[1]+0.5*k1[1])/a[1]
                       +18.0*r[2]*r[2]*alpha[2]*(Phi[2]+0.5*k1[2])/a[2]  -6.0*r[3]*r[3]*alpha[3]*(Phi[3]+0.5*k1[3])/a[3]
                            +r[4]*r[4]*alpha[4]*(Phi[4]+0.5*k1[4])/a[4])/(12.0*deltaR*r[1]*r[1]);
        for(int ir=2;ir<nR-2;ir++){
            j2[ir] = deltaT*(Pi[ir]+0.5*l1[ir])*alpha[ir]/a[ir];
            k2[ir] = deltaT*(alpha[ir-2]*(Pi[ir-2]+0.5*l1[ir-2])/a[ir-2] -8.0*alpha[ir-1]*(Pi[ir-1]+0.5*l1[ir-1])/a[ir-1]
                        +8.0*alpha[ir+1]*(Pi[ir+1]+0.5*l1[ir+1])/a[ir+1]     -alpha[ir+2]*(Pi[ir+2]+0.5*l1[ir+2])/a[ir+2])/(12.0*deltaR);
            l2[ir] = deltaT*(r[ir-2]*r[ir-2]*alpha[ir-2]*(Phi[ir-2]+0.5*k1[ir-2])/a[ir-2] -8.0*r[ir-1]*r[ir-1]*alpha[ir-1]*(Phi[ir-1]+0.5*k1[ir-1])/a[ir-1]
                        +8.0*r[ir+1]*r[ir+1]*alpha[ir+1]*(Phi[ir+1]+0.5*k1[ir+1])/a[ir+1]     -r[ir+2]*r[ir+2]*alpha[ir+2]*(Phi[ir+2]+0.5*k1[ir+2])/a[ir+2])/(12.0*deltaR*r[ir]*r[ir]);
        }
        j2[nR-2] = 0;
        k2[nR-2] = 0;
        l2[nR-2] = 0;
        j2[nR-1] = 0;
        k2[nR-1] = 0;
        l2[nR-1] = 0;

        j3[0] = deltaT*(Pi[0]+0.5*l2[0])*alpha[0]/a[0];
        k3[0] = deltaT*(-25.0*alpha[0]*(Pi[0]+0.5*l2[0])/a[0] +48.0*alpha[1]*(Pi[1]+0.5*l2[1])/a[1] 
                        -36.0*alpha[2]*(Pi[2]+0.5*l2[2])/a[2] +16.0*alpha[3]*(Pi[3]+0.5*l2[3])/a[3]
                         -3.0*alpha[4]*(Pi[4]+0.5*l2[4])/a[4])/(12.0*deltaR);
        l3[0] = deltaT*(-25.0*r[0]*r[0]*alpha[0]*(Phi[0]+0.5*k2[0])/a[0] +48.0*r[1]*r[1]*alpha[1]*(Phi[1]+0.5*k2[1])/a[1]
                        -36.0*r[2]*r[2]*alpha[2]*(Phi[2]+0.5*k2[2])/a[2] +16.0*r[3]*r[3]*alpha[3]*(Phi[3]+0.5*k2[3])/a[3]
                         -3.0*r[4]*r[4]*alpha[4]*(Phi[4]+0.5*k2[4])/a[4])/(12.0*deltaR*r[0]*r[0]);
        j3[1] = deltaT*(Pi[1]+0.5*l2[1])*alpha[1]/a[1];
        k3[1] = deltaT*(-3.0*alpha[0]*(Pi[0]+0.5*l2[0])/a[0] -10.0*alpha[1]*(Pi[1]+0.5*l2[1])/a[1]
                       +18.0*alpha[2]*(Pi[2]+0.5*l2[2])/a[2]  -6.0*alpha[3]*(Pi[3]+0.5*l2[3])/a[3]
                            +alpha[4]*(Pi[4]+0.5*l2[4])/a[4])/(12.0*deltaR);
        l3[1] = deltaT*(-3.0*r[0]*r[0]*alpha[0]*(Phi[0]+0.5*k2[0])/a[0] -10.0*r[1]*r[1]*alpha[1]*(Phi[1]+0.5*k2[1])/a[1]
                       +18.0*r[2]*r[2]*alpha[2]*(Phi[2]+0.5*k2[2])/a[2]  -6.0*r[3]*r[3]*alpha[3]*(Phi[3]+0.5*k2[3])/a[3]
                            +r[4]*r[4]*alpha[4]*(Phi[4]+0.5*k2[4])/a[4])/(12.0*deltaR*r[1]*r[1]);
        for(int ir=2;ir<nR-2;ir++){
            j3[ir] = deltaT*(Pi[ir]+0.5*l2[ir])*alpha[ir]/a[ir];
            k3[ir] = deltaT*(alpha[ir-2]*(Pi[ir-2]+0.5*l2[ir-2])/a[ir-2] -8.0*alpha[ir-1]*(Pi[ir-1]+0.5*l2[ir-1])/a[ir-1]
                        +8.0*alpha[ir+1]*(Pi[ir+1]+0.5*l2[ir+1])/a[ir+1]     -alpha[ir+2]*(Pi[ir+2]+0.5*l2[ir+2])/a[ir+2])/(12.0*deltaR);
            l3[ir] = deltaT*(r[ir-2]*r[ir-2]*alpha[ir-2]*(Phi[ir-2]+0.5*k2[ir-2])/a[ir-2] -8.0*r[ir-1]*r[ir-1]*alpha[ir-1]*(Phi[ir-1]+0.5*k2[ir-1])/a[ir-1]
                        +8.0*r[ir+1]*r[ir+1]*alpha[ir+1]*(Phi[ir+1]+0.5*k2[ir+1])/a[ir+1]     -r[ir+2]*r[ir+2]*alpha[ir+2]*(Phi[ir+2]+0.5*k2[ir+2])/a[ir+2])/(12.0*deltaR*r[ir]*r[ir]);
        }
        j3[nR-2] = 0;
        k3[nR-2] = 0;
        l3[nR-2] = 0;
        j3[nR-1] = 0;
        k3[nR-1] = 0;
        l3[nR-1] = 0;

        j4[0] = deltaT*(Pi[0]+l3[0])*alpha[0]/a[0];
        k4[0] = deltaT*(-25.0*alpha[0]*(Pi[0]+l3[0])/a[0] +48.0*alpha[1]*(Pi[1]+l3[1])/a[1] 
                        -36.0*alpha[2]*(Pi[2]+l3[2])/a[2] +16.0*alpha[3]*(Pi[3]+l3[3])/a[3]
                         -3.0*alpha[4]*(Pi[4]+l3[4])/a[4])/(12.0*deltaR);
        l4[0] = deltaT*(-25.0*r[0]*r[0]*alpha[0]*(Phi[0]+k3[0])/a[0] +48.0*r[1]*r[1]*alpha[1]*(Phi[1]+k3[1])/a[1]
                        -36.0*r[2]*r[2]*alpha[2]*(Phi[2]+k3[2])/a[2] +16.0*r[3]*r[3]*alpha[3]*(Phi[3]+k3[3])/a[3]
                         -3.0*r[4]*r[4]*alpha[4]*(Phi[4]+k3[4])/a[4])/(12.0*deltaR*r[0]*r[0]);
        j4[1] = deltaT*(Pi[1]+l3[1])*alpha[1]/a[1];
        k4[1] = deltaT*(-3.0*alpha[0]*(Pi[0]+l3[0])/a[0] -10.0*alpha[1]*(Pi[1]+l3[1])/a[1]
                       +18.0*alpha[2]*(Pi[2]+l3[2])/a[2]  -6.0*alpha[3]*(Pi[3]+l3[3])/a[3]
                            +alpha[4]*(Pi[4]+l3[4])/a[4])/(12.0*deltaR);
        l4[1] = deltaT*(-3.0*r[0]*r[0]*alpha[0]*(Phi[0]+k3[0])/a[0] -10.0*r[1]*r[1]*alpha[1]*(Phi[1]+k3[1])/a[1]
                       +18.0*r[2]*r[2]*alpha[2]*(Phi[2]+k3[2])/a[2]  -6.0*r[3]*r[3]*alpha[3]*(Phi[3]+k3[3])/a[3]
                            +r[4]*r[4]*alpha[4]*(Phi[4]+k3[4])/a[4])/(12.0*deltaR*r[1]*r[1]);
        for(int ir=2;ir<nR-2;ir++){
            j4[ir] = deltaT*(Pi[ir]+l3[ir])*alpha[ir]/a[ir];
            k4[ir] = deltaT*(alpha[ir-2]*(Pi[ir-2]+l3[ir-2])/a[ir-2] -8.0*alpha[ir-1]*(Pi[ir-1]+l3[ir-1])/a[ir-1]
                        +8.0*alpha[ir+1]*(Pi[ir+1]+l3[ir+1])/a[ir+1]     -alpha[ir+2]*(Pi[ir+2]+l3[ir+2])/a[ir+2])/(12.0*deltaR);
            l4[ir] = deltaT*(r[ir-2]*r[ir-2]*alpha[ir-2]*(Phi[ir-2]+k3[ir-2])/a[ir-2] -8.0*r[ir-1]*r[ir-1]*alpha[ir-1]*(Phi[ir-1]+k3[ir-1])/a[ir-1]
                        +8.0*r[ir+1]*r[ir+1]*alpha[ir+1]*(Phi[ir+1]+k3[ir+1])/a[ir+1]     -r[ir+2]*r[ir+2]*alpha[ir+2]*(Phi[ir+2]+k3[ir+2])/a[ir+2])/(12.0*deltaR*r[ir]*r[ir]);
        }
        j4[nR-2] = 0;
        k4[nR-2] = 0;
        l4[nR-2] = 0;
        j4[nR-1] = 0;
        k4[nR-1] = 0;
        l4[nR-1] = 0;

        for(int ir=0;ir<nR;ir++){
            phi[ir] += (j1[ir] +2.0*j2[ir] +2.0*j3[ir] +j4[ir])/6.0;
            Phi[ir] += (k1[ir] +2.0*k2[ir] +2.0*k3[ir] +k4[ir])/6.0;
             Pi[ir] += (l1[ir] +2.0*l2[ir] +2.0*l3[ir] +l4[ir])/6.0;
        }
    }

    history[0] = rHist;
    history[1] = phiHist;
    history[2] = PhiHist;
    history[3] = PiHist;
    history[4] = aHist;
    history[5] = alphaHist;
    return history;
}

void print_data(double** hist,double* model_parameters,int iterations,int maxR,double deltaR,int nP,time_t totalTime,double epsilon){
    int print_iterations, printR;
    print_iterations = iterations/SAVE_ITERATION;
    printR = (maxR/deltaR)/SAVE_RES;
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
    //Add time to filename
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char fileName[50];
    snprintf(fileName, sizeof(fileName), "Output_%02d%02d%02d.dat", tm.tm_hour, tm.tm_min, tm.tm_sec);
    FILE* data = fopen(fileName,"w");
    //Print all parameters 
    fprintf(data,"Metric: Choptuik (rewritten)\n");
    fprintf(data,"Function type: Exponential\n");
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
    double *initial_params = malloc(sizeof(float)*3);
    initial_params[0] = 0.000008;
    initial_params[1] = 20.0;
    initial_params[2] = 3.0;
    double deltaR = 0.001;
    int maxR = 100;
    int nR = maxR/deltaR;
    double **hist, **initial_func;
    double *r, *phi, *Phi, *Pi, *a, *alpha;
    int iterations = 80000;

    initial_func = initialize_field(initial_params,deltaR,maxR);
    r = initial_func[0];
    phi = initial_func[1];
    Phi = initial_func[2];
    Pi = initial_func[3];

    time_t initTime = time(NULL);
    hist = iteration(r,phi,Phi,Pi,deltaR,nR,iterations);
    time_t finalTime = time(NULL);
    time_t timeDelta = (finalTime-initTime);

    print_data(hist,initial_params,iterations,maxR,deltaR,1,timeDelta,0.0);
}