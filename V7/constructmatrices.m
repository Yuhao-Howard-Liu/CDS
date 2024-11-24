k1=[1 1 1];
sig1=[0.002 0.002 0.002];
sig2=[0.004 0.004 0.004];
theta1=[0.05 0.05 0.05];
theta2=[0.015 0.015 0.015];
kappa=[3 5];
expectedexppsi=[exp(0) exp(1*(exp(-1)-1))];
state=1;
u=[1 1 1];
alpr=[0.005 1 1 0];
alplambda=[0.005 0 1 1];
alprecover=[1 1 1 1];
alprecoverzeros=[0 0 0 0];
Q=[-0.1 0.1; 1 -1];
initialstate=1;
maturities=[1 2 3 5 7 10];

T=10;
%no of partitions for the sample Markov chain
N=101;
d=T/(N-1);
%no of Markov chains
NN=100;
%10000 MCs are good sd/mean=0.35%
%number of observation
NNN=1000;
x0=0.05*(ones(3,NNN));
%x0=0.05*rand(3,NNN);
tic
barraynorecovery=zeros(N,3,NN);
cmatrixnorecovery=zeros(N,NN);
parraynorecovery=zeros(N,3,NN);
qmatrixnorecovery=zeros(N,NN);
kappamatrixnorecovery=zeros(N,NN);
etamatrixnorecovery=zeros(N,NN);
for i=1:NN
[t,s]=ctmcgenerator(T,initialstate,Q);
    for v=1:N      

    [bbb,ccc,ppp,qqq,kkk,eee]=condexpect((v-1)*d,t,s,alpr,alplambda,alprecoverzeros,k1,sig1,sig2,theta1,theta2,kappa,expectedexppsi);

    barraynorecovery(v,:,i)=bbb;
    cmatrixnorecovery(v,i)=ccc;
    parraynorecovery(v,:,i)=ppp;
    qmatrixnorecovery(v,i)=qqq;
    kappamatrixnorecovery(v,i)=kkk;
    etamatrixnorecovery(v,i)=eee;
    end
end



tic
barraynorecovery=reshape(barraynorecovery,3,N,NN);
bx=pagemtimes(x0',barraynorecovery);
bx=reshape(bx,N,NN,NNN);
parraynorecovery=reshape(parraynorecovery,3,N,NN);
px=pagemtimes(x0',parraynorecovery);
px=reshape(px,N,NN,NNN);

pathaveraged1=reshape(mean(exp(cmatrixnorecovery+bx),2),N,NNN);
pathaveraged2=reshape(mean(kappamatrixnorecovery.*exp(cmatrixnorecovery+bx),2),N,NNN);
pathaveraged3=reshape(mean((qmatrixnorecovery+px).*kappamatrixnorecovery.*exp(cmatrixnorecovery+bx),2),N,NNN);
%payment 1.2
buyer=zeros(6,NNN);
%buyer=trapz(0:d:T,pathaveraged);
buyer(1,:)=trapz(0:d:maturities(1),pathaveraged1(1:(1+maturities(1)*1/d),:));
for i=2:6
   buyer(i,:)=buyer((i-1),:)+trapz(maturities(i-1):d:maturities(i),pathaveraged1((1+maturities(i-1)*1/d):(1+maturities(i)*1/d),:));
    
end

%payment 1.5
onepointfive=zeros(6,NNN);
onepointfive(1,:)=trapz(0:d:maturities(1),pathaveraged2(1:(1+maturities(1)*1/d),:));
for i=2:6
   onepointfive(i,:)=onepointfive((i-1),:)+trapz(maturities(i-1):d:maturities(i),pathaveraged2((1+maturities(i-1)*1/d):(1+maturities(i)*1/d),:));
    
end
onepointfive=alplambda(1).*onepointfive;

%payment 1.6
onepointsix=zeros(6,NNN);
onepointsix(1,:)=trapz(0:d:maturities(1),pathaveraged3(1:(1+maturities(1)*1/d),:));
for i=2:6
   onepointsix(i,:)=onepointsix((i-1),:)+trapz(maturities(i-1):d:maturities(i),pathaveraged3((1+maturities(i-1)*1/d):(1+maturities(i)*1/d),:));
    
end


sqrt(mean((buyer-mean(buyer,2)).^2,2))./mean(buyer,2).*100


toc