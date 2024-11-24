tic
Q=[-0.2331 0.2331; 1.0257 -1.0257].*[2 2;0.5 0.5];
initialstate=1;
sig1=[0.04 0.15 0.005];
sig2=[0.06 0.18 0.02];
theta1_p=[0.05 0.02 0.001];
theta2_p=[0.025 0.04 0.005];
k1_p=[0.1 1 1];
k2_p=[0.12 2 2];
gam1=[0.8 3 5];
gam2=[0.9 6 10];
kapparf=[2 2.5];
kappaccp=[4 5];
k1=k1_p-sig1.*gam1;
k2=k2_p-sig2.*gam2;
theta1=k1_p./k1.*theta1_p;
theta2=k2_p./k2.*theta2_p;
rtimes=1/2;
l1times=2;
l2times=5;
M=25000;
N=300;
Tmax=5;
expectedpsirf=[10/11 5/6];
expectedpsiccp=[8/9 4/5];
w20=0.3;
w21=0.4;
w10=0.24;
w11=0.22;
c=0.05;
upft=zeros(41,1);
load("xzero.mat")


for n = 1:14
x0=X(n,:);
initialstate=1;
[premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21);
upft(n,:)=protection-CCPloss-premium*c;
end

for n = 15:30
x0=X(n,:);
initialstate=2;
[premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21);
upft(n,:)=protection-CCPloss-premium*c;
end
for n = 31:41
x0=X(n,:);
initialstate=1;
[premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21);
upft(n,:)=protection-CCPloss-premium*c;
end
toc

save pathsp0.mat upft