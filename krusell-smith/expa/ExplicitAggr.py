from __future__ import division
import numpy as np
from scipy import interpolate, linalg

# Defining the transition matrix.

DurationUnempG = 1.5
DurationUnempB = 2.5
corr = .25

UnempG = .04
UnempB = .1
Unemp = np.array([[UnempG], [UnempB]])

DurationZG = 8
DurationZB = 8

pZG = 1 - 1 / DurationZG
pZB = 1 - 1 / DurationZB

PZ = np.array([[pZG, 1 - pZG], [1 - pZB, pZB]])

p22 = 1 - 1 / DurationUnempG
p21 = 1 - p22
p11 = ((1 - UnempG) - UnempG * p21) / (1 - UnempG)

P11 = np.array([[p11, 1 - p11], [p21, p22]])

p22 = 1 - 1 / DurationUnempB
p21 = 1 - p22
p11 = ((1 - UnempB) - UnempB * p21) / (1 - UnempB)

P00 = np.array([[p11, 1 - p11], [p21, p22]])

p22 = (1 + corr) * p22
p21 = 1 - p22
p11 = ((1 - UnempB) - UnempG * p21) / (1 - UnempG)

P10 = np.array([[p11, 1 - p11], [p21, p22]])

p22 = (1 - corr) * (1 - 1 / DurationUnempG)
p21 = 1 - p22
p11 = ((1 - UnempG) - UnempB * p21) / (1 - UnempB)

P01 = np.array([[p11,1 - p11], [p21, p22]])

P = np.vstack((np.hstack((PZ[0, 0] * P11, PZ[0, 1] * P10)), 
              np.hstack((PZ[1, 0] * P01, PZ[1, 1] * P00))))

# Model Parameters

alpha  = 0.36             
beta   = 0.99
delta  = .025
sigma  = 1
phi    = 0
Nk     = 250
NK     = 12
UI     = 0.15
h      = 1 / (1 - UnempB)
tau    = np.array([[UI * UnempG / (h * (1 - UnempG)), UI * UnempB / (h * (1 - UnempB))]])
BiasCorrection = 1

if BiasCorrection == 1:
    de = 0.01257504725554
    du = 0.03680683961167
else:
    de=0
    du=0

deltaZ = 0.01

# Grid

kmin   = phi
kmax   = 200

zg = 1 + deltaZ
zb = 1 - deltaZ
Z = np.array([[zg], [zb]])
ZZ = Z

Kumin = 33
Kemin = 35
Kumax = 42.5
Kemax = 43.5

kptemp = np.linspace(0, np.log(kmax + 1 - kmin), Nk).reshape((Nk, 1))
kp = np.exp(kptemp) - 1 + kmin

Ke = np.linspace(Kemin, Kemax, NK).reshape((NK, 1))
Ku = np.linspace(Kumin, Kumax, NK).reshape((NK, 1))

# Initial policy functions

# Individual policy function
kpp1 = np.empty((2, 2, NK, Nk, NK))
kpp2 = np.empty((2, 2, NK, Nk, NK))

# multi-dimensional array indexing differs between NumPy and Matlab
for ii in range(2):
    for i in range(2):
        for j in range(NK):
            for l in range(NK):
                kpp1[ii, i, l, :, j] = (1 - delta) * kp.flatten()
                kpp2[ii, i, l, :, j] = .3 * (1 - delta) * kp.flatten()

# deviates from Matlab behavior!
kmat, Kemat, Kumat = np.meshgrid(kp, Ke, Ku, indexing='xy') 

# Initial aggregate policy function (guess: unit root.)
A = np.empty((2, 2, 2))
A[0] = linalg.inv(np.array([[P11[0, 0], P11[1,0] * UnempG / (1 - UnempG)], 
                                [P11[0, 1] * (1 - UnempG) / UnempG, P11[1, 1]]]))
A[1] = linalg.inv(np.array([[P00[0,0], P00[1,0] * UnempB / (1 - UnempB)], 
                            [P00[0, 1] * (1 - UnempB)/UnempB, P00[1, 1]]]))

Kpeu = np.empty((2,1))
Kpe  = np.empty((2, NK, NK))
Kpu  = np.empty((2, NK, NK))
for ii in range(2):
    for i in range(NK):
        for j in range(NK):
            Kpeu = A[ii, :, :].dot(np.array([Ke[i], Ku[j]]))
            Kpe[ii, i, j] = Kpeu[0]
            Kpu[ii, i, j] = Kpeu[1]

KpeN = np.empty((2, 2, NK, NK))
KpeN[0, 0, :, :] = ((P11[0, 0] * (1 - UnempG) * Kpe[0, :, :] + 
                     P11[1, 0] * UnempG * Kpu[0, :, :]) / (1 - UnempG))
KpeN[1, 0, :, :] = ((P01[0, 0] * (1 - UnempB) * Kpe[1, :, :] + 
                     P01[1, 0] * UnempB * Kpu[1, :, :]) / (1 - UnempG))
KpeN[0, 1, :, :] = ((P10[0, 0] * (1 - UnempG) * Kpe[0, :, :] + 
                     P10[1, 0] * UnempG * Kpu[0, :, :]) / (1 - UnempB))
KpeN[1, 1, :, :] = ((P00[0, 0] * (1 - UnempB) * Kpe[1, :, :] + 
                     P00[1, 0] * UnempB * Kpu[1, :, :]) / (1 - UnempB))

KpuN = np.empty((2, 2, NK, NK))
KpuN[0, 0, :, :] = ((P11[0, 1] * (1 - UnempG) * Kpe[0, :, :] + 
                     P11[1, 1] * UnempG * Kpu[0, :, :]) / (UnempG))
KpuN[1, 0, :, :] = ((P01[0, 1] * (1 - UnempB) * Kpe[1, :, :] + 
                     P01[1, 1] * UnempB * Kpu[1, :, :]) / (UnempG))
KpuN[0, 1, :, :] = ((P10[0, 1] * (1 - UnempG) * Kpe[0, :, :] + 
                     P10[1, 1] * UnempG * Kpu[0, :, :]) / (UnempB))
KpuN[1, 1, :, :] = ((P00[0, 1] * (1 - UnempB) * Kpe[1, :, :] + 
                     P00[1, 1] * UnempB * Kpu[1, :, :]) / (UnempB))

KpeN = np.minimum(KpeN, np.max(Ke))
KpeN = np.maximum(KpeN, np.min(Ke))
KpuN = np.minimum(KpuN, np.max(Ku))
KpuN = np.maximum(KpuN, np.min(Ku))


# Solving the individual's Euler equation with endogenous gridpoints.

ConvCrit = 1
s = 0

print 'Solving the individual and aggregate problem until ConvCrit < 1e-6'

while s < 1:#ConvCrit > 1e-6:
    s = s + 1

    # not sure whether or not this should be inside the while loop?
    k = np.empty((2, 2, NK, Nk, NK))
    
    for i in range(2): # pages
        for j in range(NK): # rows
            for l in range(NK): # cols
                Kp = (1 - Unemp[i, 0]) * Kpe[i, j, l] + Unemp[i, 0] * Kpu[i, j, l]
                
                Rp = np.array([(1 + alpha * zg * (Kp / (h * (1 - UnempG)))**(alpha - 1) - delta), 
                               (1 + alpha * zg * (Kp / (h * (1 - UnempG)))**(alpha - 1) - delta),
                               (1 + alpha * zb * (Kp / (h * (1 - UnempB)))**(alpha - 1) - delta),
                               (1 + alpha * zb * (Kp / (h * (1 - UnempB)))**(alpha - 1) - delta)])
                
                Wp = np.array([h * (1 - alpha) * zg * (Kp / (h * (1 - UnempG)))**alpha * (1 - tau[0, 0]),
                              (1 - alpha) * zg * (Kp / (h * (1 - UnempG)))**alpha,
                              h * (1 - alpha) * zb * (Kp / (h * (1 - UnempB)))**alpha * (1 - tau[0, 1]),
                              (1 - alpha) * zb * (Kp / (h * (1 - UnempB)))**alpha])
                
                # correct dimensions are (Nk, 4)
                RHS = np.hstack((beta * Rp[0] * (Rp[0] * kp + Wp[0] - kpp1[i, 0, l, :, j].reshape((Nk, 1)))**(-sigma),
                                 beta * Rp[1] * (Rp[1] * kp + UI * Wp[1] - kpp2[i, 0, l, :, j].reshape((Nk, 1)))**(-sigma),
                                 beta * Rp[2] * (Rp[2] * kp + Wp[2] - kpp1[i, 1, l, :, j].reshape((Nk, 1)))**(-sigma),
                                 beta * Rp[3] * (Rp[3] * kp + UI * Wp[3] - kpp2[i, 1, l, :, j].reshape((Nk, 1)))**(-sigma)))
                
                C1 = (RHS.dot(P[2 * i, :].reshape((P.shape[0], 1))))**(-1 / sigma)
                C2 = (RHS.dot(P[2 * i + 1, :].reshape((P.shape[0], 1))))**(-1/sigma)
                
                K = (1 - Unemp[i, 0]) * Ke[j, 0] + Unemp[i, 0] * Ku[l, 0]
            
                k[i, 0, l, :, j] = ((C1 - h * (1 - alpha) * Z[i, 0] * (K / (h * (1 - Unemp[i, 0])))**alpha * (1 - tau[0, i]) + kp).flatten() / 
                                    ((1 + alpha * Z[i, 0] * (K / (h * (1 - Unemp[i, 0])))**(alpha - 1) - delta)))
                k[i, 1, l, :, j] = ((C2 - UI * (1 - alpha) * Z[i, 0] * (K / (h * (1 - Unemp[i, 0])))**alpha + kp).flatten() /
                                    ((1 + alpha * Z[i, 0] * (K / (h * (1 - Unemp[i, 0])))**(alpha - 1) - delta)))
    
    
    ConvCrit = 0

    kpmat = np.empty((2, 2, NK, Nk, NK))
    for i in range(2):
        for j in range(NK):
            for l in range(NK):
                x  = k[i, 0, l, :, j]
                Y  = kp.flatten()
                xi = kp.flatten()
                if np.min(k[i, 0, l, :,j]) > 0:
                    x = np.append(0, x)
                    Y = np.append(0, Y)
                    kpmat[i, 0, l, :, j] = interpolate.interp1d(x, Y, kind='linear')(xi)
                else:
                    kpmat[i, 0, l, :, j] = interpolate.interp1d(x, Y, kind='linear')(xi)
                
                x  = k[i, 1, l, :, j]
                Y  = kp.flatten()
                xi = kp.flatten()
                if np.min(k[i, 1, l, :, j]) > 0:
                    x = np.append(0, x)
                    Y  = np.append(0, Y)
                    kpmat[i, 1, l, :,j] = interpolate.interp1d(x, Y, kind='linear')(xi)
                else:
                    kpmat[i, 1, l, :, j] = interpolate.interp1d(x, Y, kind='linear')(xi)

    # Update the individual policy functions for t+1 (interpolation using griddata is unacceptably slow!):
    kpp2New = np.empty((2, 2, NK, Nk, NK))
    for ii in range(2):
        for i in range(2):
            for j in range(NK):
                for l in range(NK):
                    # interpolation points needs to be (NK * Nk * NK, 3)
                    points = np.hstack((Kemat.reshape((NK * Nk * NK, 1)),  
                                        kmat.reshape((NK * Nk * NK, 1)),
                                        Kumat.reshape((NK * Nk * NK, 1))))
                    # values needs to be (NK * Nk * NK,)
                    values = kpmat[i, 0, :, :, :].flatten()
                    # desired interpolation points needs to be (Nk, 3)
                    xi = (np.repeat(KpeN[ii, i, j, l], Nk), 
                          kp.flatten(), 
                          np.repeat(KpuN[ii, i, j, l], Nk))
                    xi = np.hstack((np.tile(KpeN[ii, i, j, l], (Nk, 1)), 
                                    kp, 
                                    np.tile(KpuN[ii, i, j, l], (Nk, 1))))
                    kpp1[ii, i, l, :, j] = interpolate.griddata(points, 
                                                                values, 
                                                                xi, 
                                                                method='linear')
                    # values needs to be (NK * Nk * NK,)
                    values = kpmat[i, 1, :, :, :].flatten()
                    kpp2New[ii, i, l, :, j] = interpolate.griddata(points, 
                                                                   values, 
                                                                   xi, 
                                                                   method='linear')
    
    # Update the convergence measure:
    ConvCrit = np.max(np.abs(kpp2New[0, 0, 0, :, 0] - kpp2[0, 0, 0, :, 0]) / (1 + np.abs(kpp2[0, 0, 0, :, 0])))
    print 'Current value of ConvCrit is', ConvCrit
     
    kpp2 = kpp2New
     
    # Update the aggregate policy function using explicit aggregation:
    KpeNew = np.empty((2, NK, NK))
    KpuNew = np.empty((2, NK, NK))
    for i in range(NK):
        for j in range(NK):
            for l in range(2):
                x = kp.flatten()
                Y = kpmat[l, 0, j, :, i]
                KpeNew[l, i, j] = interpolate.interp1d(x, Y, kind='linear')(Ke[i, 0]) + de
                Y = kpmat[l, 1, j, :, i]
                KpuNew[l, i, j] = interpolate.interp1d(x, Y, kind='linear')(Ku[j, 0]) + du

    if s > 200:
        rho = 0;
        Kpe = rho * Kpe + (1 - rho) * KpeNew
        Kpu = rho * Kpu + (1 - rho) * KpuNew
    
        KpeN[0, 0, :, :] = ((P11[0, 0] * (1 - UnempG) * Kpe[0, :, :] + 
                             P11[1, 0] * UnempG * Kpu[0, :, :]) / (1 - UnempG))
        KpeN[1, 0, :, :] = ((P01[0, 0] * (1 - UnempB) * Kpe[1, :, :] + 
                             P01[1, 0] * UnempB * Kpu[1, :, :]) / (1 - UnempG))
        KpeN[0, 1, :, :] = ((P10[0, 0] * (1 - UnempG) * Kpe[0, :, :] + 
                             P10[1, 0] * UnempG * Kpu[0, :, :]) / (1 - UnempB))
        KpeN[1, 1, :, :] = ((P00[0, 0] * (1 - UnempB) * Kpe[1, :, :] + 
                             P00[1, 0] * UnempB * Kpu[1, :, :]) / (1 - UnempB))

        KpuN[0, 0, :, :] = ((P11[0, 1] * (1 - UnempG) * Kpe[0, :, :] + 
                             P11[1, 1] * UnempG * Kpu[0, :, :]) / (UnempG))
        KpuN[1, 0, :, :] = ((P01[0, 1] * (1 - UnempB) * Kpe[1, :, :] + 
                             P01[1, 1] * UnempB * Kpu[1, :, :]) / (UnempG))
        KpuN[0, 1, :, :] = ((P10[0, 1] * (1 - UnempG) * Kpe[0, :, :] + 
                             P10[1, 1] * UnempG * Kpu[0, :, :]) / (UnempB))
        KpuN[1, 1, :, :] = ((P00[0, 1] * (1 - UnempB) * Kpe[1, :, :] + 
                             P00[1, 1] * UnempB * Kpu[1, :, :]) / (UnempB))
        
        KpeN = np.minimum(KpeN, np.max(Ke))
        KpeN = np.maximum(KpeN, np.min(Ke))
        KpuN = np.minimum(KpuN, np.max(Ku))
        KpuN = np.maximum(KpuN, np.min(Ku))

print ('The individual and aggregate problem has converged. Simulation will ' + 
       'proceed until s=10000')

# Initial Distribution

pdistyoung = np.loadtxt('pdistyoung.txt')
InDist = pdistyoung[:, 1:]

NDist = len(InDist)

kk = np.linspace(0, kmax, NDist).reshape((NDist, 1))

Pe = InDist[:, 1]
Pu = InDist[:, 0]

# Exogenous shocks

Z    = np.loadtxt('Z.txt');
ZSim = -Z + 3

KeSim = np.zeros(Z.shape)
KuSim = KeSim

KeSim[0] = (kk.T).dot(Pe)
KuSim[0] = (kk.T).dot(Pu)

KeImp = KeSim
KuImp = KuSim

KeFit = KeSim
KuFit = KuSim

KuMat2, KeMat2 = np.meshgrid(Ku, Ke)

ind_switch = np.loadtxt('ind_switch.txt')
ind_switch = -ind_switch + 3
Kind = np.zeros(Z.shape);
Cind = np.zeros((len(Z) - 1, 1))
Kind[0] = 43
percentile = np.zeros((len(Z), 6))
moments    = np.zeros((len(Z), 5))
momentse   = np.zeros((len(Z), 5))
momentsu   = np.zeros((len(Z), 5))
Rvec       = np.zeros((len(Z), 1))
Wvec       = np.zeros((len(Z), 1))

# Let's rock'n'roll

s = 0
SimLength = 10000 - 1
for i in range(SimLength):
    s = s + 1
    KpeSim = interp2(KuMat2,KeMat2,Kpe(:,:,ZSim(i)),KuSim(i),KeSim(i));
    KpuSim = interp2(KuMat2,KeMat2,Kpu(:,:,ZSim(i)),KuSim(i),KeSim(i));

    KpeFit = interp2(KuMat2,KeMat2,Kpe(:,:,ZSim(i)),KuImp(i),KeImp(i));
    KpuFit = interp2(KuMat2,KeMat2,Kpu(:,:,ZSim(i)),KuImp(i),KeImp(i));

    kprimeE = interp3(Kemat,kmat,Kumat,kpmat(:,:,:,1,ZSim(i)),kk'*Pe,kk,kk'*Pu);
    kprimeU = interp3(Kemat,kmat,Kumat,kpmat(:,:,:,2,ZSim(i)),kk'*Pe,kk,kk'*Pu);

    Kind(i+1) = interp3(Kemat,kmat,Kumat,kpmat(:,:,:,ind_switch(i),ZSim(i)),kk'*Pe,Kind(i),kk'*Pu);
    K = (1-Unemp(ZSim(i)))*KeImp(i)+Unemp(ZSim(i))*KuImp(i);
    r = 1+alpha*ZZ(ZSim(i))*(K/(h*(1-Unemp(ZSim(i)))))^(alpha-1)-delta;
    w = (1-alpha)*ZZ(ZSim(i))*(K/(h*(1-Unemp(ZSim(i)))))^(alpha);
    Cind(i,1) = r*Kind(i)+(2-ind_switch(i))*h*w*(1-tau(ZSim(i)))+(ind_switch(i)-1)*w*UI-Kind(i+1);
    Rvec(i) = r;
    Wvec(i) = w;

    kprimeE = min(kprimeE,max(kk));
    kprimeE = max(kprimeE,min(kk));

    kprimeU = min(kprimeU,max(kk));
    kprimeU = max(kprimeU,min(kk));

    P = (1-Unemp(ZSim(i)))*Pe+Unemp(ZSim(i))*Pu;
    CP = cumsum(P);
    IP5 = find(CP<0.05,1,'last');
    percentile(i,1)=(0.05-CP(IP5+1))/(CP(IP5)-CP(IP5+1))*kk(IP5)+(1-(0.05-CP(IP5+1))/(CP(IP5)-CP(IP5+1)))*kk(IP5+1);
    IP10 = find(CP<0.1,1,'last');
    percentile(i,2)=(0.1-CP(IP10+1))/(CP(IP10)-CP(IP10+1))*kk(IP10)+(1-(0.1-CP(IP10+1))/(CP(IP10)-CP(IP10+1)))*kk(IP10+1);

    CPe = cumsum(Pe);
    IPe5 = find(CPe<0.05,1,'last');
    percentile(i,3)=(0.05-CPe(IPe5+1))/(CPe(IPe5)-CPe(IPe5+1))*kk(IPe5)+(1-(0.05-CPe(IPe5+1))/(CPe(IPe5)-CPe(IPe5+1)))*kk(IPe5+1);
    IPe10 = find(CPe<0.1,1,'last');
    percentile(i,4)=(0.1-CPe(IPe10+1))/(CPe(IPe10)-CPe(IPe10+1))*kk(IPe10)+(1-(0.1-CPe(IPe10+1))/(CPe(IPe10)-CPe(IPe10+1)))*kk(IPe10+1);

    CPu = cumsum(Pu);
    IPu5 = find(CPu<0.05,1,'last');
    percentile(i,5)=(0.05-CPu(IPu5+1))/(CPu(IPu5)-CPu(IPu5+1))*kk(IPu5)+(1-(0.05-CPu(IPu5+1))/(CPu(IPu5)-CPu(IPu5+1)))*kk(IPu5+1);
    IPu10 = find(CPu<0.1,1,'last');
    percentile(i,6)=(0.1-CPu(IPu10+1))/(CPu(IPu10)-CPu(IPu10+1))*kk(IPu10)+(1-(0.1-CPu(IPu10+1))/(CPu(IPu10)-CPu(IPu10+1)))*kk(IPu10+1);

    moments(i,1) = kk'*P;
    momentse(i,1) = kk'*Pe;
    momentsu(i,1) = kk'*Pu;
    for j in range(2, 5 + 1):
        moments(i,j) = (((kk.^j)'*P))^(1/j)/moments(i,1);
        momentse(i,j) = (((kk.^j)'*Pe))^(1/j)/momentse(i,1);
        momentsu(i,j) = (((kk.^j)'*Pu))^(1/j)/momentsu(i,1);

    for j in range(NDist):
        Ie(j,1) = find(kprimeE(j)>=kk,1,'last');
        Iu(j,1) = find(kprimeU(j)>=kk,1,'last');

    Ie = min(Ie,NDist-1);
    Iu = min(Iu,NDist-1);

    rhoE = (kprimeE-kk(Ie+1))./(kk(Ie)-kk(Ie+1));

    rhoU = (kprimeU-kk(Iu+1))./(kk(Iu)-kk(Iu+1));

    Le1 = zeros(NDist,1);
    Le2 = Le1;
    Lu1 = Le1;
    Lu2 = Le1;

    for jj in range(NDist):
        Le1(Ie(jj)) = rhoE(jj)*Pe(jj)+Le1(Ie(jj));
        Le2(Ie(jj)+1) = (1-rhoE(jj))*Pe(jj)+Le2(Ie(jj)+1);
        Lu1(Iu(jj)) = rhoU(jj)*Pu(jj)+Lu1(Iu(jj));
        Lu2(Iu(jj)+1) = (1-rhoU(jj))*Pu(jj)+Lu2(Iu(jj)+1);
 
    PPe = Le1+Le2;
    PPu = Lu1+Lu2;

    if ZSim(i)==1
        if ZSim(i+1)==1
            KeSim(i+1) = (P11(1,1)*(1-UnempG)*KpeSim+P11(2,1)*UnempG*KpuSim)./(1-UnempG);
            KuSim(i+1) = (P11(1,2)*(1-UnempG)*KpeSim+P11(2,2)*UnempG*KpuSim)./(UnempG);
            KeFit(i+1) = (P11(1,1)*(1-UnempG)*KpeFit+P11(2,1)*UnempG*KpuFit)./(1-UnempG);
            KuFit(i+1) = (P11(1,2)*(1-UnempG)*KpeFit+P11(2,2)*UnempG*KpuFit)./(UnempG);
            Pe = (P11(1,1)*(1-UnempG)*PPe+P11(2,1)*UnempG*PPu)./(1-UnempG);
            Pu = (P11(1,2)*(1-UnempG)*PPe+P11(2,2)*UnempG*PPu)./(UnempG);
            KeImp(i+1) = kk'*Pe;
            KuImp(i+1) = kk'*Pu;
        end
        if ZSim(i+1)==2
            KeSim(i+1) = (P10(1,1)*(1-UnempG)*KpeSim+P10(2,1)*UnempG*KpuSim)./(1-UnempB);
            KuSim(i+1) = (P10(1,2)*(1-UnempG)*KpeSim+P10(2,2)*UnempG*KpuSim)./(UnempB);
            KeFit(i+1) = (P10(1,1)*(1-UnempG)*KpeFit+P10(2,1)*UnempG*KpuFit)./(1-UnempB);
            KuFit(i+1) = (P10(1,2)*(1-UnempG)*KpeFit+P10(2,2)*UnempG*KpuFit)./(UnempB);
            Pe = (P10(1,1)*(1-UnempG)*PPe+P10(2,1)*UnempG*PPu)./(1-UnempB);
            Pu = (P10(1,2)*(1-UnempG)*PPe+P10(2,2)*UnempG*PPu)./(UnempB);
            KeImp(i+1) = kk'*Pe;
            KuImp(i+1) = kk'*Pu;
        end
    end
    if ZSim(i)==2
        if ZSim(i+1)==1
            KeSim(i+1) = (P01(1,1)*(1-UnempB)*KpeSim+P01(2,1)*UnempB*KpuSim)./(1-UnempG);
            KuSim(i+1) = (P01(1,2)*(1-UnempB)*KpeSim+P01(2,2)*UnempB*KpuSim)./(UnempG);
            KeFit(i+1) = (P01(1,1)*(1-UnempB)*KpeFit+P01(2,1)*UnempB*KpuFit)./(1-UnempG);
            KuFit(i+1) = (P01(1,2)*(1-UnempB)*KpeFit+P01(2,2)*UnempB*KpuFit)./(UnempG);
            Pe = (P01(1,1)*(1-UnempB)*PPe+P01(2,1)*UnempB*PPu)./(1-UnempG);
            Pu = (P01(1,2)*(1-UnempB)*PPe+P01(2,2)*UnempB*PPu)./(UnempG);
            KeImp(i+1) = kk'*Pe;
            KuImp(i+1) = kk'*Pu;
        end
        if ZSim(i+1)==2
            KeSim(i+1) = (P00(1,1)*(1-UnempB)*KpeSim+P00(2,1)*UnempB*KpuSim)./(1-UnempB);
            KuSim(i+1) = (P00(1,2)*(1-UnempB)*KpeSim+P00(2,2)*UnempB*KpuSim)./(UnempB);
            KeFit(i+1) = (P00(1,1)*(1-UnempB)*KpeFit+P00(2,1)*UnempB*KpuFit)./(1-UnempB);
            KuFit(i+1) = (P00(1,2)*(1-UnempB)*KpeFit+P00(2,2)*UnempB*KpuFit)./(UnempB);
            Pe = (P00(1,1)*(1-UnempB)*PPe+P00(2,1)*UnempB*PPu)./(1-UnempB);
            Pu = (P00(1,2)*(1-UnempB)*PPe+P00(2,2)*UnempB*PPu)./(UnempB);
            KeImp(i+1) = kk'*Pe;
            KuImp(i+1) = kk'*Pu;
        end
    end
end

KSim = (1 - Unemp[ZSim]) * KeSim + Unemp[ZSim] * KuSim
KImp = (1 - Unemp[ZSim]) * KeImp + Unemp[ZSim] * KuImp
KFit = (1 - Unemp[ZSim]) * KeFit + Unemp[ZSim] * KuFit

ZZ = np.array([[zg], [zb]])
Y  = ZZ[ZSim] * (KImp**alpha) * (((1 - Unemp[ZSim]) * h)**(1 - alpha))
C  = Y[:-1] - KImp[1:] + (1 - delta) * KImp[1:-1]

Yind = Cind + Kind[1:] - (1 - delta) * Kind[:-1]

print 'Correlation of individual and aggregate consumption:', corrcoef([Cind,C])
print 'Correlation of individual consumption and aggregate income:', corrcoef([Cind,Y(1:end-1)])
print 'Correlation of individual consumption and aggregate capital:', corrcoef([Cind,KImp(1:end-1)])
print 'Correlation of individual consumption and individual income:', corrcoef([Cind,Yind])
print 'Correlation of individual consumption and individual capital:', corrcoef([Cind,Kind(1:end-1)])
print 'Standard deviation of individual consumption:', std(Cind)
print 'Standard deviation of individual capital:', std(Kind)
print 'Autocorrelation of individual consumption:', corrcoef([Cind(1:end-3),Cind(2:end-2),Cind(3:end-1),Cind(4:end)])
print 'Autocorrelation of individual capital:', corrcoef([Kind(1:end-3),Kind(2:end-2),Kind(3:end-1),Kind(4:end)])
print 'Autocorrelation of individual consumption growth:', cgrowth = log(Cind(2:end))-log(Cind(1:end-1));
corrcoef([cgrowth(1:end-3),cgrowth(2:end-2),cgrowth(3:end-1),cgrowth(4:end)])
print 'Max error Ke (%)', 100*max(abs(log(KeSim)-log(KeImp)))
print 'Max error Ku (%)', 100*max(abs(log(KuSim)-log(KuImp)))
print 'Mean error Ke (%)', 100*mean(abs(log(KeSim)-log(KeImp)))
print 'Mean error Ku (%)', 100*mean(abs(log(KuSim)-log(KuImp)))
print 'R-Square K', 1-var(KImp-KFit)/var(KImp)
print 'R-Square Ke', 1-var(KeImp-KeFit)/var(KeImp)
print 'R-Square Ku', 1-var(KuImp-KuFit)/var(KuImp)