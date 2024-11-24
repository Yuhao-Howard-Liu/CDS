tic
Q=[-0.2331 0.2331; 1.0257 -1.0257].*[2 2;0.5 0.5];
initialstate=1;
maturities=[1 5 10];
sig1=[0.04 0.15 0.01];
sig2=[0.06 0.18 0.04];
theta1_p=[0.05 0.03 0.005];
theta2_p=[0.025 0.06 0.02];
k1_p=[0.1 1 1];
k2_p=[0.12 3 3];
gam1=[0.8 4 5];
gam2=[0.9 8 10];
kapparf=[2 2.5];
kappaccp=[4 5];
k1=k1_p-sig1.*gam1;
k2=k2_p-sig2.*gam2;
theta1=k1_p./k1.*theta1_p;
theta2=k2_p./k2.*theta2_p;
rtimes=1/2;
l1times=2;
l2times=4;
x0=[0.05 0.03 0.005];
expectedpsirf=[5/6 2/3];
expectedpsiccp=[8/9 2/3];
w20=0.4;
w21=0.5;
w10=0.24;
w11=0.22;
Tmax=5;
M=1000;
N=50;
c=0.05;
upft1=zeros(41,1);
load('xzero.mat')

for n = 1:14
x0=X(n,:);
initialstate=1;
[premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21);
upft1(n,:)=CCPloss;
end

for n = 15:30
x0=X(n,:);
initialstate=2;
[premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21);
upft1(n,:)=CCPloss;
end
for n = 31:41
x0=X(n,:);
initialstate=1;
[premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21);
upft1(n,:)=CCPloss;
end
toc

save pathsp1.mat upft1