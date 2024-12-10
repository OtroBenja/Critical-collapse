#pragma once

#include "constants.h"

double get_single_mass(double *r, double *Phi, double *Pi, double *a, double maxR, double deltaR){
    double mass = 0;
    int nR = (int)(maxR/deltaR);
    for(int ir=0;ir<nR;ir++){
        mass += 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir])*deltaR;
    }
    return mass;
}

double *get_mass(double *r, double *Phi, double *Pi, double *a, double maxR, double deltaR){
    int nR = (int)(maxR/deltaR);
    double *mass = malloc(sizeof(double)*nR);
    mass[0] = 2.0*PI*(Phi[0]*Phi[0] +Pi[0]*Pi[0])*r[0]*r[0]/(a[0]*a[0])*deltaR;
    for(int ir=1;ir<nR;ir++){
        mass[ir] = mass[ir-1] + 2.0*PI*(Phi[ir]*Phi[ir] +Pi[ir]*Pi[ir])*r[ir]*r[ir]/(a[ir]*a[ir])*deltaR;
    }
    return mass;
}

#if METRIC == 0
void metric_iteration(double a0,double alpha0,double* Beta, double* Beta1_2, double* a, double* alpha, double* r, int nR, double deltaR, bool normalize){}
#endif

#if METRIC == 1
void metric_iteration(double a0,double alpha0,double* Beta, double* Beta1_2, double* a, double* alpha, double* r, int nR, double deltaR, bool normalize){
    
    double m1, n1, m2, n2, m3, n3, m4, n4, norm;
    a[0] = a0;
    alpha[0] = alpha0;
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
            a[ir+1] =     a[ir] +(m1 +2.0*(m2+m3) +m4)/6.0;
        alpha[ir+1] = alpha[ir] +(n1 +2.0*(n2+n3) +n4)/6.0;
    }
    if (normalize) {norm = 1.0/(alpha[nR-1]*a[nR-1]);
    } else {norm = 1.0;}
    for(int ir=0;ir<nR;ir++){
        alpha[ir] = alpha[ir]*norm;
    }
}
#endif