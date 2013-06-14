"""Python implementation of den Haan and Rendahl (2009) algorithm for solving 
the Krusell and Smith (1999) model.  This code is definitely not "Pythonic" and
is instead used to demonstrate how easily MatLab code can be adapted to Python.

"""
from __future__ import division
import numpy as np
from scipy.interpolate import interp1d, griddata

def interpn(*args, **kw):
    """Interpolation on N-D. 

    ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
    where the arrays x, y, z, ... define a rectangular grid
    and a.shape == (len(x), len(y), len(z), ...)
    """
    method = kw.pop('method', 'cubic')
    if kw:
        raise ValueError("Unknown arguments: " % kw.keys())
    nd = (len(args) - 1) // 2
    if len(args) != 2 * nd + 1:
        raise ValueError("Wrong number of arguments")
    q  = args[:nd]
    qi = args[nd+1:]
    a  = args[nd]
    for j in range(nd):
        a = interp1d(q[j], a, axis=j, kind=method)(qi[j])
    return a
    

########## Model Parameters ##########

# Capital's share of income
alpha = 0.36

# Capital depreciation rate 
delta = 0.025

# Individual discount factor             
beta = 0.99

# Coefficient of relative risk aversion
sigma = 1.0

# Aggregate productivity shock in the "good" state
z_G = 1.01

# Average duration of the "good" state (in quarters)
durationZG = 8.0

# Unemployment rate in the "good" state
u_G = 0.04

# Average duration of unemployment spell in "good" state (in quarters)
durationUnempG = 1.5

# Aggregate productivity shock in the "bad" state
z_B = 0.99

# Average duration of the "bad" state (in quarters)
durationZB = 8.0

# Unemployment rate in the "bad" state
u_B = 0.1

# Average duration of unemployment spell in "bad" state (in quarters)
durationUnempB = 2.5

# Correlation between aggregate shock and individual employment shocks
corr = 0.25

# fraction of wage paid when unemployed (replacement rate!)
mu = 0.15

# government taxes employed at rate tau and redistributes to unemployed
tau = np.array([(mu * u_G) * ((1 - u_B) / (1 - u_G)), mu * u_B])

########## Defining the transition matrices ##########

##### Aggregate state transition #####

# Probability that economy remains in the "good" state
probZG = 1 - 1 / durationZG

# Probability that economy remains in the "bad" state
probZB = 1 - 1 / durationZB

# Transition matrix for aggregate state of the economy
transitionZ = np.array([[probZG, 1 - probZG], 
                        [1 - probZB, probZB]])

##### Employment transitions in "good" state #####

# Probability of remaining unemployed 
probUUG = 1 - 1 / durationUnempG

# Probability of finding a job
probUEG = 1 - probUUG

# Probability of keeping a job 
probEEG = ((1 - u_G) - u_G * probUEG) / (1 - u_G)

# Probability of losing a job 
probEUG = 1 - probEEG

# Transistion matrix for employment
transitionGG = np.array([[probEEG, probEUG], 
                         [probUEG, probUUG]])

##### Employment transiton in "bad" state #####

# Probability of remaining unemployed 
probUUB = 1 - 1 / durationUnempB

# Probability of finding a job  
probUEB = 1 - probUUB

# Probability of keeping a job
probEEB = ((1 - u_B) - u_B * probUEB) / (1 - u_B)

# Probability of losing a job 
probEUB = 1 - probEEB
  
# Transition matrix for employment 
transitionBB = np.array([[probEEB, probEUB], 
                         [probUEB, probUUB]])

##### Employment transitions, "good" to "bad" #####

# Probability of remaining unemployed when switching
probUUGB = (1 + corr) * probUUB

# Probability of finding a job when switching 
probUEGB = 1 - probUUGB

# Probability of remaining employed when switching 
probEEGB = ((1 - u_B) - u_G * probUEGB) / (1 - u_G)

# Probability of losing a job when switching 
probEUGB = 1 - probEEGB

# Transition matrix for employment when switching 
transitionGB = np.array([[probEEGB, probEUGB], 
                          [probUEGB, probUUGB]])

##### Employment transitions, "bad" to "good" #####

# Probability of remaining unemployed when switching 
probUUBG = (1 - corr) * probUUG

# Probability of finding a job when switching 
probUEBG = 1 - probUUBG

# Probability of keeping a job when switching
probEEBG = ((1 - u_G) - u_B * probUEBG) / (1 - u_B)

# Probability of losing a job when switching
probEUBG = 1 - probEEBG

# Transition matrix for employment when switching
transitionBG = np.array([[probEEBG, probEUBG], 
                          [probUEBG, probUUBG]])

#####  Transition matrices for employment #####

# Note that this is a 2D array!
P = np.vstack((np.hstack((transitionZ[0,0] * transitionGG, 
                          transitionZ[0,1] * transitionGB)),
               np.hstack((transitionZ[1,0] * transitionBG, 
                          transitionZ[1,1] * transitionBB))))

########## Specifying the initial grid of state variables ##########

# Xpa implements the endogenous grid point method proposed by Carroll (2006).

##### Grid for individual capital stock k #####

# Specify the number of nodes for approximating k_prime
N_k = 250

# To get more nodes close to the constraint, use equidistant nodes in log-space
kmin   = 0.0
kmax   = 200.0

kp = kmin + np.logspace(0, np.log(1 + kmax), N_k, base=np.exp(1)) - 1 

##### Grid for aggregate capital stock of unemployed #####

# Specify the number of nodes for approximating Ku
NK     = 12

# Just use linearly spaced grid points
Kumin = 33.0
Kumax = 42.5

Ku = np.linspace(Kumin, Kumax, NK).reshape(NK, 1)

##### Grid for aggregate capital stock of employed #####

# Just use linearly spaced grid points
Kemin = 35.0
Kemax = 43.5

Ke = np.linspace(Kemin, Kemax, NK).reshape(NK, 1)

########## Initial policy functions ##########

##### Individual agent policy functions #####
""" For each set of values for employment shock, e, and aggregate productivity 
shock, z, the individual policy function is then the three dimensional 
interpolation of the optimal policy choices at the nodes of the (k, Ke, Ku) 
grid. The individual policy function is thus a five dimensional object. The 
coefficients of this optimal policy can be solved for via time iteration.

Given that employment shock can take two values, e={0, 1}, and aggregate 
productivity shock can take two values, z={z_G, z_B}, the shape of the 5D 
individual policy function should be (2, 2, NK, kp, NK).

"""

# Individual policy (conditional on being employed?)
kpp1 = np.tile((1 - delta) * kp, reps=(2, 2, NK, 1, NK))

# Individual policy (conditional on being unemployed?) Why 0.3??
kpp2 = np.tile(0.3 * (1 - delta) * kp, reps=(2, 2, NK, 1, NK))

# meshgrid3D is coming out in NumPy 1.7.1
Kemat = np.tile(Ke.T, reps=(NK, kp.shape[0], 1))
kmat  = np.tile(kp, reps=(NK, 1, Ke.shape[0]))
Kumat = np.tile(Ku[:, np.newaxis], reps=(1,kp.shape[0], 12)) 

##### Aggregate policy function (initial guess: unit root) #####

A = np.array([np.linalg.inv([[transitionGG[0,0], transitionGG[1,0] * u_G / (1 - u_G)], 
                             [transitionGG[0,1] * (1 - u_G) / u_G, transitionGG[1,1]]]),
              np.linalg.inv([[transitionBB[0,0], transitionBB[1,0] * u_B / (1 - u_B)], 
                             [transitionBB[0,1] * (1 - u_B) / u_B, transitionBB[1,1]]])])

# End of period aggregate capital stock held by employed agents
Kpe = np.zeros((2,NK,NK))

# End of period aggregate capital stock held by unemployed agents
Kpu = np.zeros((2,NK,NK))

for ii in xrange(2):
    for i in xrange(NK):
        for j in xrange(NK):
            Kpeu = A[ii,:,:].dot(np.array([Ke[i], Ku[j]]))
            Kpe[ii,i,j] = Kpeu[0]
            Kpu[ii,i,j] = Kpeu[1]

# Next period's start of period capital stock held by employed agents
KpeN = np.array([[(transitionGG[0,0] * (1 - u_G) * Kpe[0,:,:] + transitionGG[1,0] * u_G * Kpu[0,:,:]) / (1 - u_G),
                  (transitionBG[0,0] * (1 - u_B) * Kpe[1,:,:] + transitionBG[1,0] * u_B * Kpu[1,:,:]) / (1 - u_G)],
                 [(transitionGB[0,0] * (1 - u_G) * Kpe[0,:,:] + transitionGB[1,0] * u_G * Kpu[0,:,:]) / (1 - u_B),
                  (transitionBB[0,0] * (1 - u_B) * Kpe[1,:,:] + transitionBB[1,0] * u_B * Kpu[1,:,:]) / (1 - u_B)]])
# Next period's start of period capital stock held by unemployed agents
KpuN = np.array([[(transitionGG[0,1] * (1 - u_G) * Kpe[0,:,:] + transitionGG[1,1] * u_G * Kpu[0,:,:]) / u_G,
                  (transitionBG[0,1] * (1 - u_B) * Kpe[1,:,:] + transitionBG[1,1] * u_B * Kpu[1,:,:]) / u_G],
                 [(transitionGB[0,1] * (1 - u_G) * Kpe[0,:,:] + transitionGB[1,1] * u_G * Kpu[0,:,:]) / u_B,
                  (transitionBB[0,1] * (1 - u_B) * Kpe[1,:,:] + transitionBB[1,1] * u_B * Kpu[1,:,:]) / u_B]])

KpeN = np.minimum(KpeN, np.max(Ke))
KpeN = np.maximum(KpeN, np.min(Ke))
KpuN = np.minimum(KpuN, np.max(Ku))
KpuN = np.maximum(KpuN, np.min(Ku))
"""
##### Solving the individual's Euler equation with endogenous gridpoints ##### 

BiasCorrection = 1

if BiasCorrection == 1:
    de = 0.01257504725554
    du = 0.03680683961167
else:
    de = 0.0
    du = 0.0
    
ConvCrit = 1
s = 0

print 'Solving the individual and aggregate problem until ConvCrit < 1e-6'

while ConvCrit > 1e-6:
    s = s+1
    k = np.zeros((2,2,NK,kp.shape[0],NK))
    
    for i in xrange(2): # pages
        for j in xrange(NK): # rows
            for l in xrange(NK): # columns
                Kp = ((1 - Unemp[i]) * Kpe[i,j,l] + Unemp[i] * Kpu[i,j,l])
                Rp = np.hstack((1 + alpha * z_G * (Kp / (h * (1 - u_G)))**(alpha - 1) - delta,
                                1 + alpha * z_G * (Kp / (h * (1 - u_G)))**(alpha - 1) - delta,
                                1 + alpha * z_B * (Kp / (h * (1 - u_B)))**(alpha - 1) - delta,
                                1 + alpha * z_B * (Kp / (h * (1 - u_B)))**(alpha - 1) - delta))
                Wp = np.hstack((h * (1 - alpha) * z_G * (Kp / (h * (1 - u_G)))**(alpha) * (1 - tau[0]),
                                (1 - alpha) * z_G * (Kp / (h * (1 - u_G)))**(alpha),
                                h * (1 - alpha) * z_B * (Kp / (h * (1 - u_B)))**(alpha) * (1 - tau[1]),
                                (1 - alpha) * z_B * (Kp / (h * (1 - u_B)))**(alpha)))
                RHS = np.hstack((beta * Rp[0] * (Rp[0] * kp + Wp[0] - kpp1[i,0,l,:,j].reshape((kp.shape[0],1)))**(-sigma), 
                                 beta * Rp[1] * (Rp[1] * kp + UI * Wp[1] - kpp2[i,0,l,:,j].reshape((kp.shape[0],1)))**(-sigma),
                                 beta * Rp[2] * (Rp[2] * kp + Wp[2] - kpp1[i,1,l,:,j].reshape((kp.shape[0],1)))**(-sigma),
                                 beta * Rp[3] * (Rp[3] * kp + UI * Wp[3] - kpp2[i,1,l,:,j].reshape((kp.shape[0],1)))**(-sigma)))
                C1 = (RHS.dot(P[2 * i,:]))**(-1 / sigma)
                C2 = (RHS.dot(P[2 * i + 1,:]))**(-1 / sigma)
                K = (1 - Unemp[i]) * Ke[j] + Unemp[i] * Ku[l]
                k[i,0,l,:,j] = (C1 - h * (1 - alpha) * Z[i] * (K / (h * (1 - Unemp[i])))**(alpha) * (1 - tau[i]) + kp.flatten()) / \
                               ((1 + alpha * Z[i] * (K / (h * (1 - Unemp[i])))**(alpha - 1) - delta))
                k[i,1,l,:,j] = (C2 - UI * (1 - alpha) * Z[i] * (K / (h * (1 - Unemp[i])))**(alpha) + kp.flatten()) / \
                               ((1 + alpha * Z[i] * (K / (h * (1 - Unemp[i])))**(alpha - 1) - delta))

    kpmat = np.zeros((2,2,NK,kp.shape[0],NK))
    for i in xrange(2): # pages
        for j in xrange(NK): # rows
            for l in xrange(NK): # columns
                if np.min(k[i,0,l,:,j]) > 0:
                    tmp_f = interp1d(np.append(0, k[i,0,l,:,j]),
                                     np.append(0, kp), kind='linear',
                                     bounds_error=False, fill_value=np.nan)
                else:
                    tmp_f = interp1d(k[i,0,l,:,j], kp.flatten(), 
                                     kind='linear', bounds_error=False,
                                     fill_value=np.nan)
                kpmat[i,0,l,:,j] = tmp_f(kp.flatten())

                if np.min(k[i,0,l,:,j]) > 0:
                    tmp_f = interp1d(np.append(0, k[i,1,l,:,j]),
                                     np.append(0, kp), kind='linear',
                                     bounds_error=False, fill_value=np.nan)
                else:
                    tmp_f = interp1d(k[i,1,l,:,j], kp.flatten(), 
                                     kind='linear', bounds_error=False,
                                     fill_value=np.nan)
                kpmat[i,1,l,:,j] = tmp_f(kp.flatten())

    # Update the individual policy functions for t+1:
    kpp2New = np.zeros((2,2,NK,kp.shape[0],NK))
    for ii in xrange(2):
        for i in xrange(2):
            for j in xrange(NK):
                for l in xrange(NK):
                    values = kpmat[i,0,:,:,:]
                    xi = (KpeN[ii,i,l,j], kp.flatten(), KpuN[ii,i,l,j])
                    kpp1[ii,i,l,:,j] = interpn(Ke.flatten(), 
                                               kp.flatten(),
                                               Ku.flatten(),
                                               values,
                                               np.repeat(KpeN[ii,i,l,j], NK),
                                               kp.flatten(), 
                                               np.repeat(KpuN[ii,i,l,j], NK))[l,:,j]
                                     
                    values = kpmat[i,1,:,:,:]
                    kpp2New[ii,i,l,:,j] = interpn(Ke.flatten(), 
                                                  kp.flatten(),
                                                  Ku.flatten(),
                                                  values,
                                                  np.repeat(KpeN[ii,i,l,j], NK),
                                                  kp.flatten(), 
                                                  np.repeat(KpuN[ii,i,l,j], NK))[l,:,j]
                    print 'Done with l=', l
                print 'Done with j=', j
            print 'Done with i=', i
        print 'Done with ii=', ii
        
    # Update the convergence measure:
    ConvCrit = np.max(np.abs(kpp2New[0,0,0,:,0] - kpp2[0,0,0,:,0]) / (1 + np.abs(kpp2[0,0,0,:,0])))

    kpp2 = kpp2New

    # Update the aggregate policy function using explicit aggregation:
    for i in xrange(NK): # columns
        for j in xrange(NK): # pages 
            for l in xrange(2): # ??
                tmp_f = interp1d(kp, kpmat[l,0,j,:,i], kind='linear',
                                 bounds_error=False, fill_value=np.nan)
                KpeNew[l,i,j] = tmp_f(Ke[i]) + de

                tmp_f = interp1d(kp, kpmat[l,1,j,:,i], kind='linear',
                                 bounds_error=False, fill_value=np.nan)
                KpuNew[l,i,j] = tmp_f(Ku[j]) + du

    if s > 200:
        rho = 0
        Kpe = rho * Kpe + (1 - rho) * KpeNew
        Kpu = rho * Kpu + (1 - rho) * KpuNew

        KpeN = np.array([[(transitionEG[0,0] * (1 - u_G) * Kpe[0,:,:] + transitionEG[1,0] * u_G * Kpu[0,:,:]) / (1 - u_G),
                          (transitionEBG[0,0] * (1 - u_B) * Kpe[1,:,:] + transitionEBG[1,0] * u_B * Kpu[1,:,:]) / (1 - u_G)],
                         [(transitionEGB[0,0] * (1 - u_G) * Kpe[0,:,:] + transitionEGB[1,0] * u_G * Kpu[0,:,:]) / (1 - u_B),
                          (transitionEB[0,0] * (1 - u_B) * Kpe[1,:,:] + transitionEB[1,0] * u_B * Kpu[1,:,:]) / (1 - u_B)]])

        KpuN = np.array([[(transitionEG[0,1] * (1 - u_G) * Kpe[0,:,:] + transitionEG[1,1] * u_G * Kpu[0,:,:]) / u_G,
                          (transitionEBG[0,1] * (1 - u_B) * Kpe[1,:,:] + transitionEBG[1,1] * u_B * Kpu[1,:,:]) / u_G],
                         [(transitionEGB[0,1] * (1 - u_G) * Kpe[0,:,:] + transitionEGB[1,1] * u_G * Kpu[0,:,:]) / u_B,
                          (transitionEB[0,1] * (1 - u_B) * Kpe[1,:,:] + transitionEB[1,1] * u_B * Kpu[1,:,:]) / u_B]])
    
        KpeN = np.minimum(KpeN,np.max(Ke))
        KpeN = np.maximum(KpeN,np.min(Ke))
        KpuN = np.minimum(KpuN,np.max(Ku))
        KpuN = np.maximum(KpuN,min(Ku))

print 'The individual and aggregate problem has converged. Simulation will proceed until s=10000'

##### Initial Distribution #####

pdistyoung = np.loadtxt('pdistyoung.txt')
InDist = pdistyoung[:,1:]

NDist = len(InDist)

kk = np.linspace[0, kmax, NDist].reshape(NDist, 1);

Pe = InDist[:,1]
Pu = InDist[:,0]

# Exogenous shocks

load('Z.txt');
ZSim = -Z+3;

KeSim = np.zeros(np.size(Z))
KuSim = KeSim

KeSim[0] = (kk.T).dot(Pe)
KuSim[0] = (kk.T).dot(Pu)

KeImp = KeSim;
KuImp = KuSim;

KeFit = KeSim;
KuFit = KuSim;

[KuMat2,KeMat2] = meshgrid(Ku,Ke);

load('ind_switch.txt')
ind_switch = -ind_switch+3;
Kind = zeros(size(Z));
Cind = zeros(length(Z)-1,1);
Kind(1) = 43;
percentile = zeros(length(Z),6);
moments = zeros(length(Z),5);
momentse = zeros(length(Z),5);
momentsu = zeros(length(Z),5);
Rvec = zeros(length(Z),1);
Wvec = zeros(length(Z),1);

% Let's rock'n'roll

s = 0;
SimLength = 10000-1;
for i=1:SimLength
s = s+1
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
ItransitionEGB = find(CP<0.1,1,'last');
percentile(i,2)=(0.1-CP(ItransitionEGB+1))/(CP(ItransitionEGB)-CP(ItransitionEGB+1))*kk(ItransitionEGB)+(1-(0.1-CP(ItransitionEGB+1))/(CP(ItransitionEGB)-CP(ItransitionEGB+1)))*kk(ItransitionEGB+1);

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
for j=2:5
    moments(i,j) = (((kk.^j)'*P))^(1/j)/moments(i,1);
    momentse(i,j) = (((kk.^j)'*Pe))^(1/j)/momentse(i,1);
    momentsu(i,j) = (((kk.^j)'*Pu))^(1/j)/momentsu(i,1);
end

for j=1:NDist
    Ie(j,1) = find(kprimeE(j)>=kk,1,'last');
    Iu(j,1) = find(kprimeU(j)>=kk,1,'last');
end

Ie = min(Ie,NDist-1);
Iu = min(Iu,NDist-1);

rhoE = (kprimeE-kk(Ie+1))./(kk(Ie)-kk(Ie+1));

rhoU = (kprimeU-kk(Iu+1))./(kk(Iu)-kk(Iu+1));

Le1 = zeros(NDist,1);
Le2 = Le1;
Lu1 = Le1;
Lu2 = Le1;

for jj=1:NDist
    Le1(Ie(jj)) = rhoE(jj)*Pe(jj)+Le1(Ie(jj));
    Le2(Ie(jj)+1) = (1-rhoE(jj))*Pe(jj)+Le2(Ie(jj)+1);
    Lu1(Iu(jj)) = rhoU(jj)*Pu(jj)+Lu1(Iu(jj));
    Lu2(Iu(jj)+1) = (1-rhoU(jj))*Pu(jj)+Lu2(Iu(jj)+1);
end
 
PPe = Le1+Le2;
PPu = Lu1+Lu2;

if ZSim(i)==1
    if ZSim(i+1)==1
        KeSim(i+1) = (transitionEG(1,1)*(1-u_G)*KpeSim+transitionEG(2,1)*u_G*KpuSim)./(1-u_G);
        KuSim(i+1) = (transitionEG(1,2)*(1-u_G)*KpeSim+transitionEG(2,2)*u_G*KpuSim)./(u_G);
        KeFit(i+1) = (transitionEG(1,1)*(1-u_G)*KpeFit+transitionEG(2,1)*u_G*KpuFit)./(1-u_G);
        KuFit(i+1) = (transitionEG(1,2)*(1-u_G)*KpeFit+transitionEG(2,2)*u_G*KpuFit)./(u_G);
        Pe = (transitionEG(1,1)*(1-u_G)*PPe+transitionEG(2,1)*u_G*PPu)./(1-u_G);
        Pu = (transitionEG(1,2)*(1-u_G)*PPe+transitionEG(2,2)*u_G*PPu)./(u_G);
        KeImp(i+1) = kk'*Pe;
        KuImp(i+1) = kk'*Pu;
    end
    if ZSim(i+1)==2
        KeSim(i+1) = (transitionEGB(1,1)*(1-u_G)*KpeSim+transitionEGB(2,1)*u_G*KpuSim)./(1-u_B);
        KuSim(i+1) = (transitionEGB(1,2)*(1-u_G)*KpeSim+transitionEGB(2,2)*u_G*KpuSim)./(u_B);
        KeFit(i+1) = (transitionEGB(1,1)*(1-u_G)*KpeFit+transitionEGB(2,1)*u_G*KpuFit)./(1-u_B);
        KuFit(i+1) = (transitionEGB(1,2)*(1-u_G)*KpeFit+transitionEGB(2,2)*u_G*KpuFit)./(u_B);
        Pe = (transitionEGB(1,1)*(1-u_G)*PPe+transitionEGB(2,1)*u_G*PPu)./(1-u_B);
        Pu = (transitionEGB(1,2)*(1-u_G)*PPe+transitionEGB(2,2)*u_G*PPu)./(u_B);
        KeImp(i+1) = kk'*Pe;
        KuImp(i+1) = kk'*Pu;
    end
end
if ZSim(i)==2
    if ZSim(i+1)==1
        KeSim(i+1) = (transitionEBG(1,1)*(1-u_B)*KpeSim+transitionEBG(2,1)*u_B*KpuSim)./(1-u_G);
        KuSim(i+1) = (transitionEBG(1,2)*(1-u_B)*KpeSim+transitionEBG(2,2)*u_B*KpuSim)./(u_G);
        KeFit(i+1) = (transitionEBG(1,1)*(1-u_B)*KpeFit+transitionEBG(2,1)*u_B*KpuFit)./(1-u_G);
        KuFit(i+1) = (transitionEBG(1,2)*(1-u_B)*KpeFit+transitionEBG(2,2)*u_B*KpuFit)./(u_G);
        Pe = (transitionEBG(1,1)*(1-u_B)*PPe+transitionEBG(2,1)*u_B*PPu)./(1-u_G);
        Pu = (transitionEBG(1,2)*(1-u_B)*PPe+transitionEBG(2,2)*u_B*PPu)./(u_G);
        KeImp(i+1) = kk'*Pe;
        KuImp(i+1) = kk'*Pu;
    end
    if ZSim(i+1)==2
        KeSim(i+1) = (transitionEB(1,1)*(1-u_B)*KpeSim+transitionEB(2,1)*u_B*KpuSim)./(1-u_B);
        KuSim(i+1) = (transitionEB(1,2)*(1-u_B)*KpeSim+transitionEB(2,2)*u_B*KpuSim)./(u_B);
        KeFit(i+1) = (transitionEB(1,1)*(1-u_B)*KpeFit+transitionEB(2,1)*u_B*KpuFit)./(1-u_B);
        KuFit(i+1) = (transitionEB(1,2)*(1-u_B)*KpeFit+transitionEB(2,2)*u_B*KpuFit)./(u_B);
        Pe = (transitionEB(1,1)*(1-u_B)*PPe+transitionEB(2,1)*u_B*PPu)./(1-u_B);
        Pu = (transitionEB(1,2)*(1-u_B)*PPe+transitionEB(2,2)*u_B*PPu)./(u_B);
        KeImp(i+1) = kk'*Pe;
        KuImp(i+1) = kk'*Pu;
    end
end


end

KSim = (1-Unemp(ZSim)).*KeSim+Unemp(ZSim).*KuSim;
KImp = (1-Unemp(ZSim)).*KeImp+Unemp(ZSim).*KuImp;
KFit = (1-Unemp(ZSim)).*KeFit+Unemp(ZSim).*KuFit;

ZZ = [z_G;z_B];
Y = ZZ(ZSim).*(KImp.^alpha).*(((1-Unemp(ZSim))*h).^(1-alpha));
C = Y(1:end-1)-KImp(2:end)+(1-delta)*KImp(1:end-1);

Yind = Cind+Kind(2:end)-(1-delta)*Kind(1:end-1);

disp('Correlation of individual and aggregate consumption:')
corrcoef([Cind,C])
disp('Correlation of individual consumption and aggregate income:')
corrcoef([Cind,Y(1:end-1)])
disp('Correlation of individual consumption and aggregate capital:')
corrcoef([Cind,KImp(1:end-1)])
disp('Correlation of individual consumption and individual income:')
corrcoef([Cind,Yind])
disp('Correlation of individual consumption and individual capital:')
corrcoef([Cind,Kind(1:end-1)])
disp('Standard deviation of individual consumption:')
std(Cind)
disp('Standard deviation of individual capital:')
std(Kind)
disp('Autocorrelation of individual consumption:')
corrcoef([Cind(1:end-3),Cind(2:end-2),Cind(3:end-1),Cind(4:end)])
disp('Autocorrelation of individual capital:')
corrcoef([Kind(1:end-3),Kind(2:end-2),Kind(3:end-1),Kind(4:end)])
disp('Autocorrelation of individual consumption growth:')
cgrowth = log(Cind(2:end))-log(Cind(1:end-1));
corrcoef([cgrowth(1:end-3),cgrowth(2:end-2),cgrowth(3:end-1),cgrowth(4:end)])
disp('Max error Ke (%)')
100*max(abs(log(KeSim)-log(KeImp)))
disp('Max error Ku (%)')
100*max(abs(log(KuSim)-log(KuImp)))
disp('Mean error Ke (%)')
100*mean(abs(log(KeSim)-log(KeImp)))
disp('Mean error Ku (%)')
100*mean(abs(log(KuSim)-log(KuImp)))
disp('R-Square K')
1-var(KImp-KFit)/var(KImp)
disp('R-Square Ke')
1-var(KeImp-KeFit)/var(KeImp)
disp('R-Square Ku')
1-var(KuImp-KuFit)/var(KuImp)
"""
