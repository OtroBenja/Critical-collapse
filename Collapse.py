import numpy as np
from numba import njit
from scipy.integrate import RK45, DOP853

deltaR = 0.01
deltaT = deltaR/5.0
iterations = 7810
min_iterations = 0

def initialize_r(r0,maxR):
    R = np.arange(0,maxR,deltaR)
    R[0] = r0
    return R

def initialize_fields(R, p0 = 0.000008):
    phi = np.empty_like(R)
    Phi = np.empty_like(R)
    Pi  = np.zeros_like(R)

    phi = p0*(R**3)*(np.e**(-((R-20)/3)**2))

    #Calcular derivada
    Phi[   0] = 0
    Phi[   1] = ( phi[  4] -6.0*phi[   3] +18.0*phi[   2] -10.0*phi[ 1] -3.0*phi[ 0])/(12.0*deltaR)
    Phi[2:-2] = ( phi[:-4] -8.0*phi[1:-3]  +8.0*phi[3:-1]      -phi[4:]             )/(12.0*deltaR)
    Phi[  -2] = (-phi[ -5] +6.0*phi[  -4] -18.0*phi[  -3] +10.0*phi[-2] +3.0*phi[-1])/(12.0*deltaR)
    Phi[  -1] = 0

    return phi, Phi, Pi

@njit
def iterate_metric(R,Phi,Pi,a,alpha):
    for nR in range(len(R)-1):
        m1 = deltaR*(    a[nR]*(2*np.pi*R[nR]*(Phi[nR]**2 + Pi[nR]**2) - 0.5*(a[nR]**2 - 1)/R[nR]))
        n1 = deltaR*(alpha[nR]*(2*np.pi*R[nR]*(Phi[nR]**2 + Pi[nR]**2) + 0.5*(a[nR]**2 - 1)/R[nR]))
        m2 = deltaR*(    a[nR]+0.5*m1)*(2*np.pi*(R[nR]+0.5*deltaR)*((0.5*Phi[nR]+0.5*Phi[nR+1])**2 + (0.5*Pi[nR]+0.5*Pi[nR+1])**2) - 0.5*((a[nR]+0.5*m1)**2 - 1)/(R[nR]+0.5*deltaR))
        n2 = deltaR*(alpha[nR]+0.5*n1)*(2*np.pi*(R[nR]+0.5*deltaR)*((0.5*Phi[nR]+0.5*Phi[nR+1])**2 + (0.5*Pi[nR]+0.5*Pi[nR+1])**2) + 0.5*((a[nR]+0.5*m1)**2 - 1)/(R[nR]+0.5*deltaR))
        m3 = deltaR*(    a[nR]+0.5*m2)*(2*np.pi*(R[nR]+0.5*deltaR)*((0.5*Phi[nR]+0.5*Phi[nR+1])**2 + (0.5*Pi[nR]+0.5*Pi[nR+1])**2) - 0.5*((a[nR]+0.5*m2)**2 - 1)/(R[nR]+0.5*deltaR))
        n3 = deltaR*(alpha[nR]+0.5*n2)*(2*np.pi*(R[nR]+0.5*deltaR)*((0.5*Phi[nR]+0.5*Phi[nR+1])**2 + (0.5*Pi[nR]+0.5*Pi[nR+1])**2) + 0.5*((a[nR]+0.5*m2)**2 - 1)/(R[nR]+0.5*deltaR))
        m4 = deltaR*(    a[nR]+m3)*(2*np.pi*(R[nR]+deltaR)*((Phi[nR+1])**2 + (Pi[nR+1])**2) - 0.5*((a[nR]+m2)**2 - 1)/(R[nR]+deltaR))
        n4 = deltaR*(alpha[nR]+n3)*(2*np.pi*(R[nR]+deltaR)*((Phi[nR+1])**2 + (Pi[nR+1])**2) + 0.5*((a[nR]+m2)**2 - 1)/(R[nR]+deltaR))

        a[nR+1]     = a[nR]     + (m1+2.0*m2+2.0*m3+m4)/6.0
        alpha[nR+1] = alpha[nR] + (n1+2.0*n2+2.0*n3+n4)/6.0

@njit
def iterate(R,Phi,Pi,Phi_hist,A_hist,iter,min_iter):
    a = np.ones_like(R)
    alpha = np.ones_like(R)
    k1 = np.empty_like(R)
    l1 = np.empty_like(R)
    k2 = np.empty_like(R)
    l2 = np.empty_like(R)
    k3 = np.empty_like(R)
    l3 = np.empty_like(R)
    k4 = np.empty_like(R)
    l4 = np.empty_like(R)

    for i in range(iter):

        # a and alpha with RK4
        iterate_metric(R,Phi,Pi,a,alpha)

        if(i>=min_iter):
            Phi_hist[i-min_iter] = Phi
            A_hist[i-min_iter] = a

        # next Phi and Pi step with RK4
        k1[0] = deltaT*(-25.0*alpha[0]*Pi[0]/a[0] +48.0*alpha[1]*Pi[1]/a[1] 
                        -36.0*alpha[2]*Pi[2]/a[2] +16.0*alpha[3]*Pi[3]/a[3]
                         -3.0*alpha[4]*Pi[4]/a[4])/(12.0*deltaR)
        l1[0] = deltaT*(-25.0*R[0]**2*alpha[0]*Phi[0]/a[0] +48.0*R[1]**2*alpha[1]*Phi[1]/a[1]
                        -36.0*R[2]**2*alpha[2]*Phi[2]/a[2] +16.0*R[3]**2*alpha[3]*Phi[3]/a[3]
                         -3.0*R[4]**2*alpha[4]*Phi[4]/a[4])/(12.0*deltaR*R[0]**2)
        k1[1] = deltaT*(-3.0*alpha[0]*Pi[0]/a[0] -10.0*alpha[1]*Pi[1]/a[1]
                       +18.0*alpha[2]*Pi[2]/a[2]  -6.0*alpha[3]*Pi[3]/a[3]
                            +alpha[4]*Pi[4]/a[4])/(12.0*deltaR)
        l1[1] = deltaT*(-3.0*R[0]**2*alpha[0]*Pi[0]/a[0] -10.0*R[1]**2*alpha[1]*Pi[1]/a[1]
                       +18.0*R[2]**2*alpha[2]*Pi[2]/a[2]  -6.0*R[3]**2*alpha[3]*Pi[3]/a[3]
                            +R[4]**2*alpha[4]*Pi[4]/a[4])/(12.0*deltaR*R[1]**2)
        k1[2:-2] = deltaT*(alpha[ :-4]*Pi[ :-4]/a[ :-4] -8.0*alpha[1:-3]*Pi[1:-3]/a[1:-3]
                      +8.0*alpha[3:-1]*Pi[3:-1]/a[3:-1]     -alpha[4:  ]*Pi[4:  ]/a[4:  ])/(12.0*deltaR)
        l1[2:-2] = deltaT*(R[ :-4]**2*alpha[ :-4]*Phi[ :-4]/a[ :-4] -8.0*R[1:-3]**2*alpha[1:-3]*Phi[1:-3]/a[1:-3]
                      +8.0*R[3:-1]**2*alpha[3:-1]*Phi[3:-1]/a[3:-1]     -R[4:  ]**2*alpha[4:  ]*Phi[4:  ]/a[4:  ])/(12.0*deltaR*R[2:-2]**2)
        k1[-2] = 0
        l1[-1] = 0

        k2[0] = deltaT*(-25.0*alpha[0]*(Pi[0]+0.5*l1[0])/a[0] +48.0*alpha[1]*(Pi[1]+0.5*l1[1])/a[1]
                        -36.0*alpha[2]*(Pi[2]+0.5*l1[2])/a[2] +16.0*alpha[3]*(Pi[3]+0.5*l1[3])/a[3]
                         -3.0*alpha[4]*(Pi[4]+0.5*l1[4])/a[4])/(12.0*deltaR)
        l2[0] = deltaT*(-25.0*R[0]**2*alpha[0]*(Phi[0]+0.5*k1[0])/a[0] +48.0*R[1]**2*alpha[1]*(Phi[1]+0.5*k1[1])/a[1]
                        -36.0*R[2]**2*alpha[2]*(Phi[2]+0.5*k1[2])/a[2] +16.0*R[3]**2*alpha[3]*(Phi[3]+0.5*k1[3])/a[3]
                         -3.0*R[4]**2*alpha[4]*(Phi[4]+0.5*k1[4])/a[4])/(12.0*deltaR*R[0]**2)
        k2[1] = deltaT*( -3.0*alpha[0]*(Pi[0]+0.5*l1[0])/a[0] -10.0*alpha[1]*(Pi[1]+0.5*l1[1])/a[1]
                        +18.0*alpha[2]*(Pi[2]+0.5*l1[2])/a[2]  -6.0*alpha[3]*(Pi[3]+0.5*l1[3])/a[3]
                             +alpha[4]*(Pi[4]+0.5*l1[4])/a[4])/(12.0*deltaR)
        l2[1] = deltaT*( -3.0*R[0]**2*alpha[0]*(Phi[0]+0.5*k1[0])/a[0] -10.0*R[1]**2*alpha[1]*(Phi[1]+0.5*k1[1])/a[1]
                        +18.0*R[2]**2*alpha[2]*(Phi[2]+0.5*k1[2])/a[2]  -6.0*R[3]**2*alpha[3]*(Phi[3]+0.5*k1[3])/a[3]
                             +R[4]**2*alpha[4]*(Phi[4]+0.5*k1[4])/a[4])/(12.0*deltaR*R[1]**2)
        k2[2:-2] = deltaT*(alpha[ :-4]*(Pi[ :-4]+0.5*l1[ :-4])/a[ :-4] -8.0*alpha[1:-3]*(Pi[1:-3]+0.5*l1[1:-3])/a[1:-3]
                      +8.0*alpha[3:-1]*(Pi[3:-1]+0.5*l1[3:-1])/a[3:-1]     -alpha[4:  ]*(Pi[4:  ]+0.5*l1[4:  ])/a[4:  ])/(12.0*deltaR)
        l2[2:-2] = deltaT*(R[ :-4]**2*alpha[ :-4]*(Phi[ :-4]+0.5*k1[ :-4])/a[ :-4] -8.0*R[1:-3]**2*alpha[1:-3]*(Phi[1:-3]+0.5*k1[1:-3])/a[1:-3]
                      +8.0*R[3:-1]**2*alpha[3:-1]*(Phi[3:-1]+0.5*k1[3:-1])/a[3:-1]     -R[4:  ]**2*alpha[4:  ]*(Phi[4:  ]+0.5*k1[4:  ])/a[4:  ])/(12.0*deltaR*R[2:-2]**2)
        k2[-2] = 0
        l2[-1] = 0

        k3[0] = deltaT*(-25.0*alpha[0]*(Pi[0]+0.5*l2[0])/a[0] +48.0*alpha[1]*(Pi[1]+0.5*l2[1])/a[1]
                        -36.0*alpha[2]*(Pi[2]+0.5*l2[2])/a[2] +16.0*alpha[3]*(Pi[3]+0.5*l2[3])/a[3]
                         -3.0*alpha[4]*(Pi[4]+0.5*l2[4])/a[4])/(12.0*deltaR)
        l3[0] = deltaT*(-25.0*R[0]**2*alpha[0]*(Phi[0]+0.5*k2[0])/a[0] +48.0*R[1]**2*alpha[1]*(Phi[1]+0.5*k2[1])/a[1]
                        -36.0*R[2]**2*alpha[2]*(Phi[2]+0.5*k2[2])/a[2] +16.0*R[3]**2*alpha[3]*(Phi[3]+0.5*k2[3])/a[3]
                         -3.0*R[4]**2*alpha[4]*(Phi[4]+0.5*k2[4])/a[4])/(12.0*deltaR*R[0]**2)
        k3[1] = deltaT*( -3.0*alpha[0]*(Pi[0]+0.5*l2[0])/a[0] -10.0*alpha[1]*(Pi[1]+0.5*l2[1])/a[1]
                        +18.0*alpha[2]*(Pi[2]+0.5*l2[2])/a[2]  -6.0*alpha[3]*(Pi[3]+0.5*l2[3])/a[3]
                             +alpha[4]*(Pi[4]+0.5*l2[4])/a[4])/(12.0*deltaR)
        l3[1] = deltaT*( -3.0*R[0]**2*alpha[0]*(Phi[0]+0.5*k2[0])/a[0] -10.0*R[1]**2*alpha[1]*(Phi[1]+0.5*k2[1])/a[1]
                        +18.0*R[2]**2*alpha[2]*(Phi[2]+0.5*k2[2])/a[2]  -6.0*R[3]**2*alpha[3]*(Phi[3]+0.5*k2[3])/a[3]
                             +R[4]**2*alpha[4]*(Phi[4]+0.5*k2[4])/a[4])/(12.0*deltaR*R[1]**2)
        k3[2:-2] = deltaT*(alpha[ :-4]*(Pi[ :-4]+0.5*l2[ :-4])/a[ :-4] -8.0*alpha[1:-3]*(Pi[1:-3]+0.5*l2[1:-3])/a[1:-3]
                      +8.0*alpha[3:-1]*(Pi[3:-1]+0.5*l2[3:-1])/a[3:-1]     -alpha[4:  ]*(Pi[4:  ]+0.5*l2[4:  ])/a[4:  ])/(12.0*deltaR)
        l3[2:-2] = deltaT*(R[ :-4]**2*alpha[ :-4]*(Phi[ :-4]+0.5*k2[ :-4])/a[ :-4] -8.0*R[1:-3]**2*alpha[1:-3]*(Phi[1:-3]+0.5*k2[1:-3])/a[1:-3]
                      +8.0*R[3:-1]**2*alpha[3:-1]*(Phi[3:-1]+0.5*k2[3:-1])/a[3:-1]     -R[4:  ]**2*alpha[4:  ]*(Phi[4:  ]+0.5*k2[4:  ])/a[4:  ])/(12.0*deltaR*R[2:-2]**2)
        k3[-2] = 0
        l3[-1] = 0

        k4[0] = deltaT*(-25.0*alpha[0]*(Pi[0]+l3[0])/a[0] +48.0*alpha[1]*(Pi[1]+l3[1])/a[1]
                        -36.0*alpha[2]*(Pi[2]+l3[2])/a[2] +16.0*alpha[3]*(Pi[3]+l3[3])/a[3]
                         -3.0*alpha[4]*(Pi[4]+l3[4])/a[4])/(12.0*deltaR)
        l4[0] = deltaT*(-25.0*R[0]**2*alpha[0]*(Phi[0]+k3[0])/a[0] +48.0*R[1]**2*alpha[1]*(Phi[1]+k3[1])/a[1]
                        -36.0*R[2]**2*alpha[2]*(Phi[2]+k3[2])/a[2] +16.0*R[3]**2*alpha[3]*(Phi[3]+k3[3])/a[3]
                         -3.0*R[4]**2*alpha[4]*(Phi[4]+k3[4])/a[4])/(12.0*deltaR*R[0]**2)
        k4[1] = deltaT*( -3.0*alpha[0]*(Pi[0]+l3[0])/a[0] -10.0*alpha[1]*(Pi[1]+l3[1])/a[1]
                        +18.0*alpha[2]*(Pi[2]+l3[2])/a[2]  -6.0*alpha[3]*(Pi[3]+l3[3])/a[3]
                             +alpha[4]*(Pi[4]+l3[4])/a[4])/(12.0*deltaR)
        l4[1] = deltaT*( -3.0*R[0]**2*alpha[0]*(Phi[0]+k3[0])/a[0] -10.0*R[1]**2*alpha[1]*(Phi[1]+k3[1])/a[1]
                        +18.0*R[2]**2*alpha[2]*(Phi[2]+k3[2])/a[2]  -6.0*R[3]**2*alpha[3]*(Phi[3]+k3[3])/a[3]
                             +R[4]**2*alpha[4]*(Phi[4]+k3[4])/a[4])/(12.0*deltaR*R[1]**2)
        k4[2:-2] = deltaT*(alpha[ :-4]*(Pi[ :-4]+l3[ :-4])/a[ :-4] -8.0*alpha[1:-3]*(Pi[1:-3]+l3[1:-3])/a[1:-3]
                      +8.0*alpha[3:-1]*(Pi[3:-1]+l3[3:-1])/a[3:-1]     -alpha[4:  ]*(Pi[4:  ]+l3[4:  ])/a[4:  ])/(12.0*deltaR)
        l4[2:-2] = deltaT*(R[ :-4]**2*alpha[ :-4]*(Phi[ :-4]+k3[ :-4])/a[ :-4] -8.0*R[1:-3]**2*alpha[1:-3]*(Phi[1:-3]+k3[1:-3])/a[1:-3]
                      +8.0*R[3:-1]**2*alpha[3:-1]*(Phi[3:-1]+k3[3:-1])/a[3:-1]     -R[4:  ]**2*alpha[4:  ]*(Phi[4:  ]+k3[4:  ])/a[4:  ])/(12.0*deltaR*R[2:-2]**2)
        k4[-2] = 0
        l4[-1] = 0

        Phi = Phi + (k1+2.0*k2+2.0*k3+k4)/6.0
        Pi  = Pi + (l1+2.0*l2+2.0*l3+l4)/6.0

def iterate_scipy(R,Phi,Pi,Phi_hist,A_hist,iterations,min_iter):
    a = np.ones_like(R)
    alpha = np.ones_like(R)
    sizeR = len(R)
    maxT = deltaT*iterations

    #define vectorized sistem of equations F
    y = np.array((Phi,Pi))
    f = np.empty_like(y)
    def F(t,Y):
        iterate_metric(R,Y[0],Y[1],a,alpha)

        f[0][0] = (-25.0*alpha[0]*Y[1][0]/a[0] +48.0*alpha[1]*Y[1][1]/a[1] 
                    -36.0*alpha[2]*Y[1][2]/a[2] +16.0*alpha[3]*Y[1][3]/a[3]
                    -3.0*alpha[4]*Y[1][4]/a[4])/(12.0*deltaR)
        f[1][0] = (-25.0*R[0]**2*alpha[0]*Y[0][0]/a[0] +48.0*R[1]**2*alpha[1]*Y[0][1]/a[1]
                    -36.0*R[2]**2*alpha[2]*Y[0][2]/a[2] +16.0*R[3]**2*alpha[3]*Y[0][3]/a[3]
                    -3.0*R[4]**2*alpha[4]*Y[0][4]/a[4])/(12.0*deltaR*R[0]**2)
        f[0][1] = (-3.0*alpha[0]*Y[1][0]/a[0] -10.0*alpha[1]*Y[1][1]/a[1]
                    +18.0*alpha[2]*Y[1][2]/a[2]  -6.0*alpha[3]*Y[1][3]/a[3]
                        +alpha[4]*Y[1][4]/a[4])/(12.0*deltaR)
        f[1][1] = (-3.0*R[0]**2*alpha[0]*Y[0][0]/a[0] -10.0*R[1]**2*alpha[1]*Y[0][1]/a[1]
                    +18.0*R[2]**2*alpha[2]*Y[0][2]/a[2]  -6.0*R[3]**2*alpha[3]*Y[0][3]/a[3]
                        +R[4]**2*alpha[4]*Y[0][4]/a[4])/(12.0*deltaR*R[1]**2)
        f[0][2:-2] = (alpha[ :-4]*Y[1][ :-4]/a[ :-4] -8.0*alpha[1:-3]*Y[1][1:-3]/a[1:-3]
                    +8.0*alpha[3:-1]*Y[1][3:-1]/a[3:-1]     -alpha[4:  ]*Y[1][4:  ]/a[4:  ])/(12.0*deltaR)
        f[1][2:-2] = (R[ :-4]**2*alpha[ :-4]*Y[0][ :-4]/a[ :-4] -8.0*R[1:-3]**2*alpha[1:-3]*Y[0][1:-3]/a[1:-3]
                    +8.0*R[3:-1]**2*alpha[3:-1]*Y[0][3:-1]/a[3:-1]     -R[4:  ]**2*alpha[4:  ]*Y[0][4:  ]/a[4:  ])/(12.0*deltaR*R[2:-2]**2)
        f[0][-2] = 0
        f[1][-2] = 0
        f[0][-1] = 0
        f[1][-1] = 0
        return f

    #if(i>=min_iter):
    #    Phi_hist[i-min_iter] = Phi
    #    A_hist[i-min_iter] = a

    # next Phi and Pi step with RK4
    RK45(F,t0=0,y0=y,t_bound=iterations,max_step=deltaT,vectorized=False,first_step=deltaT)
        
    

r = initialize_r(deltaR*0.1,maxR = 100)
phi, Phi, Pi = initialize_fields(r)

Phi_hist = np.empty((iterations-min_iterations,len(r)))
A_hist   = np.empty((iterations-min_iterations,len(r)))

print('Compiling functions...')
iterate_scipy(r,Phi,Pi,Phi_hist,A_hist,2,1)
print('Iteration started')
iterate_scipy(r,Phi,Pi,Phi_hist,A_hist,iterations,min_iterations)
print('Iteration finished')

#Write data to files
Phi_file = open('Phi.bin','wb')
A_file = open('A.bin','wb')
R_file = open('R.bin','wb')
Phi_hist = Phi_hist.flatten()
A_hist = A_hist.flatten()
Phi_file.write(bytes(Phi_hist))
A_file.write(bytes(A_hist))
R_file.write(bytes(r))
Phi_file.close()
A_file.close()
R_file.close()
print('Data written to files')
