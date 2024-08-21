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
#endif

#define MINKOWSKI 0
#define SAVE_RES 500
#define SAVE_ITERATION 1000
#define ITERATIONS 20000
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
    for(int i=0;i<nR;i++)
        r[i] = i*deltaR;
    if(fType==0){
        for(int i=0;i<nR;i++)
            phi[i] = p0*tanh((i*deltaR-r0)/d);
    }
    if(fType==1){
        for(int i=0;i<nR;i++)
            phi[i] = p0*pow(i*deltaR,3)*pow(E,-pow((i*deltaR-r0)/d,q));
    }

    //Calculate initial Phi
    Phi[0] = 0;
    Phi[nR-1] = 0;
    for(int i=1;i<nR-1;i++){
        Phi[i] = 0.5*(phi[i+1]-phi[i-1])/deltaR;
    }

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

double** minkowski_iteration(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration){

    int nR = maxR/deltaR;
    double deltaT=deltaR/5.;
    double* k1 = malloc(sizeof(double)*nR);
    double* l1 = malloc(sizeof(double)*nR);
    double* k2 = malloc(sizeof(double)*nR);
    double* l2 = malloc(sizeof(double)*nR);
    double* k3 = malloc(sizeof(double)*nR);
    double* l3 = malloc(sizeof(double)*nR);
    double* k4 = malloc(sizeof(double)*nR);
    double* l4 = malloc(sizeof(double)*nR);
    double* Xhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Yhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Fhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Rhistory = malloc(sizeof(double)*(nR/SAVE_RES));
    double** hist = malloc(sizeof(double*)*4);
    int save_count = save_iteration;

    for(int i=0;i<iterations;i++){
        //Save values of X and Y
        if(save_count == save_iteration){
            //printf("iteration %d\n",i);
            #pragma omp parallel for
            for(int ir=0;ir<(nR/SAVE_RES);ir++){
                //Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]*Phi[ir*SAVE_RES]*sqrt(2*PI)/a;
                //Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]* Pi[ir*SAVE_RES]*sqrt(2*PI)/a;
                Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = Phi[ir*SAVE_RES];
                Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =  Pi[ir*SAVE_RES];
                Fhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = phi[ir*SAVE_RES];
            }
            save_count=0;
        }
        save_count+=1;
        //Advance Pi and Phi using RK4
        //calculate k1 and l1
        k1[0] = 0.5*(-3*Pi[0] +4*Pi[1] -Pi[2])/deltaT;
        k1[nR-1] = 0;
        l1[0] = 0.5*(-3*r[0]*r[0]*Phi[0] +4*r[1]*r[1]*Phi[1] -r[2]*r[2]*Phi[2])/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l1[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k1[ir] = 0.5*(Pi[ir+1] -Pi[ir-1])/deltaT;
            l1[ir] = 0.5*(r[ir+1]*r[ir+1]*Phi[ir+1] -r[ir-1]*r[ir-1]*Phi[ir-1])/(deltaT*r[ir]*r[ir]);
        }
        //calculate k2 and l2
        k2[0] = 0.5*(-3*(Pi[0]+0.5*l1[0]*deltaT) +4*(Pi[1]+0.5*l1[1]*deltaT) -(Pi[2]+0.5*l1[2]*deltaT))/deltaT;
        k2[nR-1] = 0;
        l2[0] = 0.5*(-3*r[0]*r[0]*(Phi[0]+0.5*k1[0]*deltaT) +4*r[1]*r[1]*(Phi[1]+0.5*k1[1]*deltaT) -r[2]*r[2]*(Phi[2]+0.5*k1[2]*deltaT))/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l2[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k2[ir] = 0.5*((Pi[ir+1]+0.5*l1[ir+1]*deltaT) -(Pi[ir-1]+0.5*l1[ir-1]*deltaT))/deltaT;
            l2[ir] = 0.5*(r[ir+1]*r[ir+1]*(Phi[ir+1]+0.5*k1[ir+1]*deltaT) -r[ir-1]*r[ir-1]*(Phi[ir-1]+0.5*k1[ir-1]*deltaT))/(deltaT*r[ir]*r[ir]);
        }
        //calculate k3 and l3
        k3[0] = 0.5*(-3*(Pi[0]+0.5*l2[0]*deltaT) +4*(Pi[1]+0.5*l2[1]*deltaT) -(Pi[2]+0.5*l2[2]*deltaT))/deltaT;
        k3[nR-1] = 0;
        l3[0] = 0.5*(-3*r[0]*r[0]*(Phi[0]+0.5*k2[0]*deltaT) +4*r[1]*r[1]*(Phi[1]+0.5*k2[1]*deltaT) -r[2]*r[2]*(Phi[2]+0.5*k2[2]*deltaT))/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l3[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k3[ir] = 0.5*((Pi[ir+1]+0.5*l2[ir+1]*deltaT) -(Pi[ir-1]+0.5*l2[ir-1]*deltaT))/deltaT;
            l3[ir] = 0.5*(r[ir+1]*r[ir+1]*(Phi[ir+1]+0.5*k2[ir+1]*deltaT) -r[ir-1]*r[ir-1]*(Phi[ir-1]+0.5*k2[ir-1]*deltaT))/(deltaT*r[ir]*r[ir]);
        }
        //calculate k4 and l4
        k4[0] = 0.5*(-3*(Pi[0]+l3[0]*deltaT) +4*(Pi[1]+l3[1]*deltaT) -(Pi[2]+l3[2]*deltaT))/deltaT;
        k4[nR-1] = 0;
        l4[0] = 0.5*(-3*r[0]*r[0]*(Phi[0]+k3[0]*deltaT) +4*r[1]*r[1]*(Phi[1]+k3[1]*deltaT) -r[2]*r[2]*(Phi[2]+k3[2]*deltaT))/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l4[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k4[ir] = 0.5*((Pi[ir+1]+l3[ir+1]*deltaT) -(Pi[ir-1]+l3[ir-1]*deltaT))/deltaT;
            l4[ir] = 0.5*(r[ir+1]*r[ir+1]*(Phi[ir+1]+k3[ir+1]*deltaT) -r[ir-1]*r[ir-1]*(Phi[ir-1]+k3[ir-1]*deltaT))/(deltaT*r[ir]*r[ir]);
        }
        //Calculate phi, Phi and Pi on next step
        #pragma omp parallel for
        for(int ir=0;ir<nR;ir++){
            Phi[ir] += (k1[ir]+2*k2[ir]+2*k3[ir]+k4[ir])*deltaT/6.;
            Pi[ir]  += (l1[ir]+2*l2[ir]+2*l3[ir]+l4[ir])*deltaT/6.;
        }
        for(int ir=nR-2;ir>-1;ir+=-1){
            phi[ir] = phi[ir+1] - Phi[ir]*deltaR; //phi is only calculated for visualization purposes
        }
    }
    for(int ir=0;ir<(nR/SAVE_RES);ir++){
        Rhistory[ir] = r[ir*SAVE_RES];
    }
    hist[0] = Rhistory;
    hist[1] = Fhistory;
    hist[2] = Xhistory;
    hist[3] = Yhistory;
    return hist;
}

double** iteration(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration){

    int nR = maxR/deltaR;
    double deltaT = deltaR/5.;
    double* a = malloc(sizeof(double)*nR);
    double* alpha = malloc(sizeof(double)*nR);
    double* Beta = malloc(sizeof(double)*nR);
    double* Beta1_2 = malloc(sizeof(double)*nR);
    double* Gamma = malloc(sizeof(double)*nR);
    double* Epsilon = malloc(sizeof(double)*nR);
    double* r2 = malloc(sizeof(double)*nR);
    double* k1 = malloc(sizeof(double)*nR);
    double* l1 = malloc(sizeof(double)*nR);
    double* k2 = malloc(sizeof(double)*nR);
    double* l2 = malloc(sizeof(double)*nR);
    double* k3 = malloc(sizeof(double)*nR);
    double* l3 = malloc(sizeof(double)*nR);
    double* k4 = malloc(sizeof(double)*nR);
    double* l4 = malloc(sizeof(double)*nR);
    double m1;
    double n1;
    double m2;
    double n2;
    double m3;
    double n3;
    double m4;
    double n4;
    double* Xhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Yhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Fhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Rhistory = malloc(sizeof(double)*(nR/SAVE_RES));
    double** hist = malloc(sizeof(double*)*4);
    int save_count = save_iteration;

    a[0] = 1.;
    alpha[0] = 1.;

    for(int ir=0;ir<nR;ir++){
        a[ir] = 1.;
        alpha[ir] = 1.;
    }

    for(int i=0;i<iterations;i++){
        
        //Solve a and alpha with known values of PHI and PI using RK4
        //Useful auxiliar variable
        #pragma omp parallel for
        for(int ir=0;ir<nR;ir++)
            Beta[ir] = 2.0*PI*r[ir]*(Pi[ir]*Pi[ir]+Phi[ir]*Phi[ir]);
        #pragma omp parallel for
        for(int ir=0;ir<nR-1;ir++)
            Beta1_2[ir] = 0.5*PI*(r[ir]+deltaR/2)*((Pi[ir]+Pi[ir+1])*(Pi[ir]+Pi[ir+1])+(Phi[ir]+Phi[ir+1])*(Phi[ir]+Phi[ir+1]));
        
        for(int ir=0;ir<nR-1;ir++){
            //calculate m1 and n1
            m1 = a[ir]*(Beta[ir]-0.5*(a[ir]*a[ir]-1)/(r[ir]+deltaR/10));
            n1 = alpha[ir]*(Beta[ir]+0.5*(a[ir]*a[ir]-1)/(r[ir]+deltaR/10));
            //calculate m2 and n2
            m2 = (a[ir]+deltaR*m1/2)*(Beta1_2[ir]-0.5*((a[ir]+deltaR*m1/2)*(a[ir]+deltaR*m1/2)-1)/(r[ir]+deltaR/2));
            n2 = (alpha[ir]+deltaR*n1/2)*(Beta1_2[ir]+0.5*((a[ir]+deltaR*m1/2)*(a[ir]+deltaR*m1/2)-1)/(r[ir]+deltaR/2));
            //calculate m3 and n3
            m3 = (a[ir]+deltaR*m2/2)*(Beta1_2[ir]-0.5*((a[ir]+deltaR*m2/2)*(a[ir]+deltaR*m2/2)-1)/(r[ir]+deltaR/2));
            n3 = (alpha[ir]+deltaR*n2/2)*(Beta1_2[ir]+0.5*((a[ir]+deltaR*m2/2)*(a[ir]+deltaR*m2/2)-1)/(r[ir]+deltaR/2));
            //calculate m4 and n4
            m4 = (a[ir]+deltaR*m3)*(Beta[ir+1]-0.5*((a[ir]+deltaR*m3)*(a[ir]+deltaR*m3)-1)/r[ir+1]);
            n4 = (alpha[ir]+deltaR*n3)*(Beta[ir+1]+0.5*((a[ir]+deltaR*m3)*(a[ir]+deltaR*m3)-1)/r[ir+1]);
            //Calculate next step for a and alpha
            a[ir+1] = a[ir] + deltaR*(m1 +2*m2 +2*m3 +m4)/6;
            alpha[ir+1] = alpha[ir] + deltaR*(n1 +2*n2 +2*n3 +n4)/6;
        }

        //Save values of X and Y
        if(save_count == save_iteration){
            //printf("iteration %d\n",i);
            #pragma omp parallel for
            for(int ir=0;ir<(nR/SAVE_RES);ir++){
                //Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]*Phi[ir*SAVE_RES]*sqrt(2*PI)/a[ir*SAVE_RES];
                //Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]* Pi[ir*SAVE_RES]*sqrt(2*PI)/a[ir*SAVE_RES];
                Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = Phi[ir*SAVE_RES];
                Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =  Pi[ir*SAVE_RES];
                Fhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = phi[ir*SAVE_RES];
                //Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = alpha[ir*SAVE_RES];
                //Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = a[ir*SAVE_RES];
            }
            save_count=0;
        }
        save_count+=1;
        //Advance Pi and Phi using RK4 
        //Usefull auxiliar variables
        #pragma omp parallel for
        for(int ir=0;ir<nR;ir++){
            r2[ir] = r[ir]*r[ir];
            Gamma[ir] = alpha[ir]/a[ir];
            Epsilon[ir] = r[ir]*r[ir]*alpha[ir]/a[ir];
        }
        //calculate k1 and l1
        k1[0] = 0.5*(-3*Gamma[0]*Pi[0] +4*Gamma[1]*Pi[1] -Gamma[2]*Pi[2])/deltaT;
        //k1[nR-1] = -0.5*(3*Phi[nR-1] -4*Phi[nR-2] +Phi[nR-3])/deltaT;
        k1[nR-1] = 0;
        l1[0] = 0.5*(-3*Epsilon[0]*Phi[0] +4*Epsilon[1]*Phi[1] -Epsilon[2]*Phi[2])/(deltaT*(r2[0]+0.01*deltaR));
        //l1[nR-1] = -0.5*(3*Epsilon[nR-1]*Gamma[nR-1]*Pi[nR-1] -4*Epsilon[nR-2]*Gamma[nR-2]*Pi[nR-2] +Epsilon[nR-3]*Gamma[nR-3]*Pi[nR-3])/(deltaT*r2[nR-1]);
        l1[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k1[ir] = 0.5*(Gamma[ir+1]*Pi[ir+1] -Gamma[ir-1]*Pi[ir-1])/deltaT;
            l1[ir] = 0.5*(Epsilon[ir+1]*Phi[ir+1] -Epsilon[ir-1]*Phi[ir-1])/(deltaT*r2[ir]);
        }
        //calculate k2 and l2
        k2[0] = 0.5*(-3*Gamma[0]*(Pi[0]+0.5*l1[0]*deltaT) +4*Gamma[1]*(Pi[1]+0.5*l1[1]*deltaT) -Gamma[2]*(Pi[2]+0.5*l1[2]*deltaT))/deltaT;
        //k2[nR-1] = -0.5*(3*(Phi[nR-1]+0.5*k1[nR-1]*deltaT) -4*(Phi[nR-2]+0.5*k1[nR-2]*deltaT) +(Phi[nR-3]+0.5*k1[nR-3]*deltaT))/deltaT;
        k2[nR-1] = 0;
        l2[0] = 0.5*(-3*Epsilon[0]*(Phi[0]+0.5*k1[0]*deltaT) +4*Epsilon[1]*(Phi[1]+0.5*k1[1]*deltaT) -Epsilon[2]*(Phi[2]+0.5*k1[2]*deltaT))/(deltaT*(r2[0]+0.01*deltaR));
        //l2[nR-1] = -0.5*(3*Epsilon[nR-1]*Gamma[nR-1]*(Pi[nR-1]+0.5*l1[nR-1]*deltaT) -4*Epsilon[nR-2]*Gamma[nR-2]*(Pi[nR-2]+0.5*l1[nR-2]*deltaT) +Epsilon[nR-3]*Gamma[nR-3]*(Pi[nR-3]+0.5*l1[nR-3]*deltaT))/(deltaT*r2[nR-1]);
        l2[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k2[ir] = 0.5*(Gamma[ir+1]*(Pi[ir+1]+0.5*l1[ir+1]*deltaT) -Gamma[ir-1]*(Pi[ir-1]+0.5*l1[ir-1]*deltaT))/deltaT;
            l2[ir] = 0.5*(Epsilon[ir+1]*(Phi[ir+1]+0.5*k1[ir+1]*deltaT) -Epsilon[ir-1]*(Phi[ir-1]+0.5*k1[ir-1]*deltaT))/(deltaT*r2[ir]);
        }
        //calculate k3 and l3
        k3[0] = 0.5*(-3*Gamma[0]*(Pi[0]+0.5*l2[0]*deltaT) +4*Gamma[1]*(Pi[1]+0.5*l2[1]*deltaT) -Gamma[2]*(Pi[2]+0.5*l2[2]*deltaT))/deltaT;
        //k3[nR-1] = -0.5*(3*(Phi[nR-1]+0.5*k2[nR-1]*deltaT) -4*(Phi[nR-2]+0.5*k2[nR-2]*deltaT) +(Phi[nR-3]+0.5*k2[nR-3]*deltaT))/deltaT;
        k3[nR-1] = 0;
        l3[0] = 0.5*(-3*Epsilon[0]*(Phi[0]+0.5*k2[0]*deltaT) +4*Epsilon[1]*(Phi[1]+0.5*k2[1]*deltaT) -Epsilon[2]*(Phi[2]+0.5*k2[2]*deltaT))/(deltaT*(r2[0]+0.01*deltaR));
        //l3[nR-1] = -0.5*(3*Epsilon[nR-1]*Gamma[nR-1]*(Pi[nR-1]+0.5*l2[nR-1]*deltaT) -4*Epsilon[nR-2]*Gamma[nR-2]*(Pi[nR-2]+0.5*l2[nR-2]*deltaT) +Epsilon[nR-3]*Gamma[nR-3]*(Pi[nR-3]+0.5*l2[nR-3]*deltaT))/(deltaT*r2[nR-1]);
        l3[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k3[ir] = 0.5*(Gamma[ir+1]*(Pi[ir+1]+0.5*l2[ir+1]*deltaT) -Gamma[ir-1]*(Pi[ir-1]+0.5*l2[ir-1]*deltaT))/deltaT;
            l3[ir] = 0.5*(Epsilon[ir+1]*(Phi[ir+1]+0.5*k2[ir+1]*deltaT) -Epsilon[ir-1]*(Phi[ir-1]+0.5*k2[ir-1]*deltaT))/(deltaT*r2[ir]);
        }
        //calculate k4 and l4
        k4[0] = 0.5*(-3*Gamma[0]*(Pi[0]+l3[0]*deltaT) +4*Gamma[1]*(Pi[1]+l3[1]*deltaT) -Gamma[2]*(Pi[2]+l3[2]*deltaT))/deltaT;
        //k4[nR-1] = -0.5*(3*(Phi[nR-1]+k3[nR-1]*deltaT) -4*(Phi[nR-2]+k3[nR-2]*deltaT) +(Phi[nR-3]+k3[nR-3]*deltaT))/deltaT;
        k4[nR-1] = 0;
        l4[0] = 0.5*(-3*Epsilon[0]*(Phi[0]+k3[0]*deltaT) +4*Epsilon[1]*(Phi[1]+k3[1]*deltaT) -Epsilon[2]*(Phi[2]+k3[2]*deltaT))/(deltaT*(r[0]+0.01*deltaR));
        //l4[nR-1] = -0.5*(3*Epsilon[nR-1]*Gamma[nR-1]*(Pi[nR-1]+l3[nR-1]*deltaT) -4*Epsilon[nR-2]*Gamma[nR-2]*(Pi[nR-2]+l3[nR-2]*deltaT) +Epsilon[nR-3]*Gamma[nR-3]*(Pi[nR-3]+l3[nR-3]*deltaT))/(deltaT*r2[nR-1]);
        l4[nR-1] = 0;
        #pragma omp parallel for
        for(int ir=1;ir<nR-1;ir++){
            k4[ir] = 0.5*(Gamma[ir+1]*(Pi[ir+1]+l3[ir+1]*deltaT) -Gamma[ir-1]*(Pi[ir-1]+l3[ir-1]*deltaT))/deltaT;
            l4[ir] = 0.5*(Epsilon[ir+1]*(Phi[ir+1]+k3[ir+1]*deltaT) -Epsilon[ir-1]*(Phi[ir-1]+k3[ir-1]*deltaT))/(deltaT*r2[ir]);
        }
        //Calculate phi, Phi and Pi on next step
        #pragma omp parallel for
        for(int ir=0;ir<nR;ir++){
            Phi[ir] += (k1[ir]+2*k2[ir]+2*k3[ir]+k4[ir])*deltaT/6.;
            Pi[ir] += (l1[ir]+2*l2[ir]+2*l3[ir]+l4[ir])*deltaT/6.;
        }
        for(int ir=nR-2;ir>-1;ir+=-1){
            phi[ir] = phi[ir+1] - Phi[ir]*deltaR; //phi is only calculated for visualization purposes
        }
    }
    for(int ir=0;ir<(nR/SAVE_RES);ir++){
        Rhistory[ir] = r[ir*SAVE_RES];
    }
    hist[0] = Rhistory;
    hist[1] = Fhistory;
    hist[2] = Xhistory;
    hist[3] = Yhistory;
    return hist;
}


void print_data(double** hist,int fType,double* model_parameters,int iterations,int maxR,double deltaR,int nP,time_t totalTime){
    int print_iterations = iterations/SAVE_ITERATION;
    int printR = (maxR/deltaR)/SAVE_RES;
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
    double q  = model_parameters[3];
    //Add time to filename
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char fileName[50];
    snprintf(fileName, sizeof(fileName), "Output_%02d%02d%02d.dat", tm.tm_hour, tm.tm_min, tm.tm_sec);
    FILE* data = fopen(fileName,"w");

    //Print all parameters
    if(MINKOWSKI) fprintf(data,"Metric: Minkowski\n");
    else          fprintf(data,"Metric: Choptuik\n" );
    if(fType) fprintf(data,"Function type: Exponential\n");
    else fprintf(data,"Function type: Hyperbolic Tan\n");
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
    for(int i=0;i<print_iterations-1;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[3][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[3][i*printR+printR-1]);
    }
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(data,"%lf,",hist[3][(print_iterations-1)*printR+ir]);
    }
    fprintf(data,"%lf",hist[3][print_iterations*printR-1]);
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
    int maxR = 70;
    if((argc>6) && atoi(argv[6])) maxR = atoi(argv[6]);
    int iterations = ITERATIONS;
    if((argc>7) && atoi(argv[7])) iterations = atoi(argv[7]);
    printf("Total iterations: %d\n",iterations);

    initial_conditions = initialize_field(fType,model_parameters,deltaR,maxR);
    r = initial_conditions[0];
    phi = initial_conditions[1];
    Phi = initial_conditions[2];
    Pi = initial_conditions[3];

    //Pass initial conditions to iteration
    if((argc>8) && atoi(argv[8])) omp_set_num_threads(atoi(argv[8]));
    time_t initTime = time(NULL);
    if(MINKOWSKI) hist = minkowski_iteration(r,phi,Phi,Pi,deltaR,maxR,iterations,SAVE_ITERATION);
    else          hist =           iteration(r,phi,Phi,Pi,deltaR,maxR,iterations,SAVE_ITERATION);
    time_t finalTime = time(NULL);
    int nP = omp_get_max_threads();
    float timeDelta = (finalTime-initTime);

    //Print simulation history to a file
    print_data(hist,fType,model_parameters,iterations,maxR,deltaR,nP,timeDelta);
    printf("Finished\n");
}




