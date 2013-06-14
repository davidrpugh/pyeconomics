clear;

% Defining the transition matrix.

DurationUnempG = 1.5;
DurationUnempB = 2.5;
corr = .25;

UnempG = .04;
UnempB = .1;
Unemp = [UnempG; UnempB];

DurationZG = 8;
DurationZB = 8;

pZG = 1-1/DurationZG;
pZB = 1-1/DurationZB;

PZ = [pZG, 1-pZG; 1-pZB, pZB];

p22 = 1-1/DurationUnempG;
p21 = 1-p22;
p11 = ((1-UnempG)-UnempG*p21)/(1-UnempG);

P11 = [p11,1-p11;p21,p22];

p22 = 1-1/DurationUnempB;
p21 = 1-p22;
p11 = ((1-UnempB)-UnempB*p21)/(1-UnempB);

P00 = [p11,1-p11;p21,p22];

p22 = (1+corr)*p22;
p21 = 1-p22;
p11 = ((1-UnempB)-UnempG*p21)/(1-UnempG);

P10 = [p11,1-p11;p21,p22];

p22 = (1-corr)*(1-1/DurationUnempG);
p21 = 1-p22;
p11 = ((1-UnempG)-UnempB*p21)/(1-UnempB);

P01 = [p11,1-p11;p21,p22];

P = [PZ(1,1)*P11,PZ(1,2)*P10;PZ(2,1)*P01,PZ(2,2)*P00];

% Model Parameters

alpha  = 0.36;             
beta   = 0.99;
delta  = .025;
sigma  = 1;
phi    = 0;
Nk     = 250;
NK     = 12;
UI     = 0.15;
h      = 1/(1-UnempB);
tau    = [UI*UnempG/(h*(1-UnempG)),UI*UnempB/(h*(1-UnempB))];
BiasCorrection = 1;

if BiasCorrection==1
    de = 0.01257504725554;
    du = 0.03680683961167;
else
    de=0;
    du=0;
end

deltaZ = 0.01;

% Grid

kmin   = phi;
kmax   = 200;

zg = 1+deltaZ;
zb = 1-deltaZ;
Z = [zg;zb];
ZZ = Z;

Kumin = 33;
Kemin = 35;
Kumax = 42.5;
Kemax = 43.5;

kptemp = linspace(0,log(kmax+1-kmin),Nk)';
kp = exp(kptemp)-1+kmin;

Ke = linspace(Kemin,Kemax,NK)';
Ku = linspace(Kumin,Kumax,NK)';

% Initial policy functions

% Individual policy function
for ii=1:2
    for i=1:2
        for j=1:NK
            for l=1:NK
                kpp1(:,j,l,i,ii) = (1-delta)*kp;
                kpp2(:,j,l,i,ii) = .3*(1-delta)*kp;
            end
        end
    end
end

[Kemat,kmat,Kumat] = meshgrid(Ke,kp,Ku);

% Initial aggregate policy function (guess: unit root.)

A(:,:,1) = inv([P11(1,1), P11(2,1)*UnempG/(1-UnempG); P11(1,2)*(1-UnempG)/UnempG, P11(2,2)]);
A(:,:,2) = inv([P00(1,1), P00(2,1)*UnempB/(1-UnempB); P00(1,2)*(1-UnempB)/UnempB, P00(2,2)]);

for ii=1:2
    for i=1:NK
        for j=1:NK
            Kpeu = A(:,:,ii)*[Ke(i);Ku(j)];
            Kpe(i,j,ii) = Kpeu(1);
            Kpu(i,j,ii) = Kpeu(2);
        end
    end
end

KpeN(:,:,1,1) = (P11(1,1)*(1-UnempG)*Kpe(:,:,1)+P11(2,1)*UnempG*Kpu(:,:,1))./(1-UnempG);
KpeN(:,:,1,2) = (P01(1,1)*(1-UnempB)*Kpe(:,:,2)+P01(2,1)*UnempB*Kpu(:,:,2))./(1-UnempG);
KpeN(:,:,2,1) = (P10(1,1)*(1-UnempG)*Kpe(:,:,1)+P10(2,1)*UnempG*Kpu(:,:,1))./(1-UnempB);
KpeN(:,:,2,2) = (P00(1,1)*(1-UnempB)*Kpe(:,:,2)+P00(2,1)*UnempB*Kpu(:,:,2))./(1-UnempB);

KpuN(:,:,1,1) = (P11(1,2)*(1-UnempG)*Kpe(:,:,1)+P11(2,2)*UnempG*Kpu(:,:,1))./(UnempG);
KpuN(:,:,1,2) = (P01(1,2)*(1-UnempB)*Kpe(:,:,2)+P01(2,2)*UnempB*Kpu(:,:,2))./(UnempG);
KpuN(:,:,2,1) = (P10(1,2)*(1-UnempG)*Kpe(:,:,1)+P10(2,2)*UnempG*Kpu(:,:,1))./(UnempB);
KpuN(:,:,2,2) = (P00(1,2)*(1-UnempB)*Kpe(:,:,2)+P00(2,2)*UnempB*Kpu(:,:,2))./(UnempB);

KpeN = min(KpeN,max(Ke));
KpeN = max(KpeN,min(Ke));
KpuN = min(KpuN,max(Ku));
KpuN = max(KpuN,min(Ku));

% Solving the individual's Euler equation with endogenous gridpoints.

ConvCrit = 1;
s = 0;

disp('Solving the individual and aggregate problem until ConvCrit<1e-6')

while ConvCrit>1e-6

s = s+1;

for i=1:2
    for j=1:NK
        for l=1:NK
            Kp = ((1-Unemp(i))*Kpe(j,l,i)+Unemp(i)*Kpu(j,l,i));
            Rp = [(1+alpha*zg*(Kp/(h*(1-UnempG)))^(alpha-1)-delta),(1+alpha*zg*(Kp/(h*(1-UnempG)))^(alpha-1)-delta), ...
                (1+alpha*zb*(Kp/(h*(1-UnempB)))^(alpha-1)-delta),(1+alpha*zb*(Kp/(h*(1-UnempB)))^(alpha-1)-delta)];
            Wp = [h*(1-alpha)*zg*(Kp/(h*(1-UnempG)))^(alpha)*(1-tau(1)),(1-alpha)*zg*(Kp/(h*(1-UnempG)))^(alpha), ...
                h*(1-alpha)*zb*(Kp/(h*(1-UnempB)))^(alpha)*(1-tau(2)),(1-alpha)*zb*(Kp/(h*(1-UnempB)))^(alpha)];
            RHS = [beta*Rp(1)*(Rp(1)*kp+Wp(1)-kpp1(:,j,l,1,i)).^(-sigma), ...
                beta*Rp(2)*(Rp(2)*kp+UI*Wp(2)-kpp2(:,j,l,1,i)).^(-sigma), ...
                beta*Rp(3)*(Rp(3)*kp+Wp(3)-kpp1(:,j,l,2,i)).^(-sigma), ...
                beta*Rp(4)*(Rp(4)*kp+UI*Wp(4)-kpp2(:,j,l,2,i)).^(-sigma)];
            C1 = (RHS*P(i*2-1,:)').^(-1/sigma);
            C2 = (RHS*P(i*2,:)').^(-1/sigma);
            K = (1-Unemp(i))*Ke(j)+Unemp(i)*Ku(l);
            k(:,j,l,1,i) = (C1-h*(1-alpha)*Z(i)*(K/(h*(1-Unemp(i))))^(alpha)*(1-tau(i))+kp) ...
                ./((1+alpha*Z(i)*(K/(h*(1-Unemp(i))))^(alpha-1)-delta));
            k(:,j,l,2,i) = (C2-UI*(1-alpha)*Z(i)*(K/(h*(1-Unemp(i))))^(alpha)+kp) ...
                ./((1+alpha*Z(i)*(K/(h*(1-Unemp(i))))^(alpha-1)-delta));
        end
    end
end
ConvCrit = 0

%{
for i=1:2
    for j=1:NK
        for l=1:NK
            if min(k(:,j,l,1,i))>0
                kpmat(:,j,l,1,i) = interp1([0;k(:,j,l,1,i)],[0;kp],kp,[],'extrap');
            else
                kpmat(:,j,l,1,i) = interp1(k(:,j,l,1,i),kp,kp,[],'extrap');
            end
            if min(k(:,j,l,2,i))>0
                kpmat(:,j,l,2,i) = interp1([0;k(:,j,l,2,i)],[0;kp],kp,[],'extrap');
            else
                kpmat(:,j,l,2,i) = interp1(k(:,j,l,2,i),kp,kp,[],'extrap');
            end
        end
    end
end

% Update the individual policy functions for t+1:
for ii=1:2
    for i=1:2
        for j=1:NK
            for l=1:NK
                kpp1(:,j,l,i,ii) = interp3(Kemat,kmat,Kumat,kpmat(:,:,:,1,i),KpeN(j,l,i,ii),kp,KpuN(j,l,i,ii));
                kpp2New(:,j,l,i,ii) = interp3(Kemat,kmat,Kumat,kpmat(:,:,:,2,i),KpeN(j,l,i,ii),kp,KpuN(j,l,i,ii));
            end
        end
    end
end

% Update the convergence measure:
ConvCrit = max(abs(kpp2New(:,1,1,1,1)-kpp2(:,1,1,1,1))./(1+abs(kpp2(:,1,1,1,1))))

kpp2 = kpp2New;

% Update the aggregate policy function using explicit aggregation:
for i=1:NK
    for j=1:NK
        for l=1:2
            KpeNew(i,j,l) = interp1(kp,kpmat(:,i,j,1,l),Ke(i))+de;
            KpuNew(i,j,l) = interp1(kp,kpmat(:,i,j,2,l),Ku(j))+du;
        end
    end
end

if s>200
    rho = 0;
    Kpe = rho*Kpe+(1-rho)*KpeNew;
    Kpu = rho*Kpu+(1-rho)*KpuNew;

    KpeN(:,:,1,1) = (P11(1,1)*(1-UnempG)*Kpe(:,:,1)+P11(2,1)*UnempG*Kpu(:,:,1))./(1-UnempG);
    KpeN(:,:,1,2) = (P01(1,1)*(1-UnempB)*Kpe(:,:,2)+P01(2,1)*UnempB*Kpu(:,:,2))./(1-UnempG);
    KpeN(:,:,2,1) = (P10(1,1)*(1-UnempG)*Kpe(:,:,1)+P10(2,1)*UnempG*Kpu(:,:,1))./(1-UnempB);
    KpeN(:,:,2,2) = (P00(1,1)*(1-UnempB)*Kpe(:,:,2)+P00(2,1)*UnempB*Kpu(:,:,2))./(1-UnempB);

    KpuN(:,:,1,1) = (P11(1,2)*(1-UnempG)*Kpe(:,:,1)+P11(2,2)*UnempG*Kpu(:,:,1))./(UnempG);
    KpuN(:,:,1,2) = (P01(1,2)*(1-UnempB)*Kpe(:,:,2)+P01(2,2)*UnempB*Kpu(:,:,2))./(UnempG);
    KpuN(:,:,2,1) = (P10(1,2)*(1-UnempG)*Kpe(:,:,1)+P10(2,2)*UnempG*Kpu(:,:,1))./(UnempB);
    KpuN(:,:,2,2) = (P00(1,2)*(1-UnempB)*Kpe(:,:,2)+P00(2,2)*UnempB*Kpu(:,:,2))./(UnempB);
    
    KpeN = min(KpeN,max(Ke));
    KpeN = max(KpeN,min(Ke));
    KpuN = min(KpuN,max(Ku));
    KpuN = max(KpuN,min(Ku));
end

% And repeat this procedure until convergence.
end
disp('The individual and aggregate problem has converged. Simulation will proceed until s=10000')

% Initial Distribution

load('pdistyoung.txt')
InDist = pdistyoung(:,2:3);

NDist = length(InDist);

kk = linspace(0,kmax,NDist)';

Pe = InDist(:,2);
Pu = InDist(:,1);

% Exogenous shocks

load('Z.txt');
ZSim = -Z+3;

KeSim = zeros(size(Z));
KuSim = KeSim;

KeSim(1) = kk'*Pe;
KuSim(1) = kk'*Pu;

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

KSim = (1-Unemp(ZSim)).*KeSim+Unemp(ZSim).*KuSim;
KImp = (1-Unemp(ZSim)).*KeImp+Unemp(ZSim).*KuImp;
KFit = (1-Unemp(ZSim)).*KeFit+Unemp(ZSim).*KuFit;

ZZ = [zg;zb];
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
%}
