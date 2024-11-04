#pragma once

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