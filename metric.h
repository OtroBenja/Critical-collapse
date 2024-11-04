
void metric_iteration(double* Beta, double* Beta1_2, double* a, double* alpha, double* r, int nR, double deltaR){
    
    double m1, n1, m2, n2, m3, n3, m4, n4, norm;
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
            a[ir+1] =     a[ir] +(m1 +2.0*(m2+m3) +m4)/6.0;
        alpha[ir+1] = alpha[ir] +(n1 +2.0*(n2+n3) +n4)/6.0;
    }
    norm = 1.0/(alpha[nR-1]*a[nR-1]);
    for(int ir=0;ir<nR;ir++){
        alpha[ir] = alpha[ir]*norm;
    }
}