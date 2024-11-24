function [cond2a, cond2b] = condexpect2(v,t,s,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsiccp,w20,w21)
%v scalar -- the discretized time points
%t vector -- the times for the MC relization
%s vector -- the states for the MC relization
%cond2a is for (6.2) in v38
%cond2b is for (6.3) in v38
t=t(t<=v);
s=s(t<=v);

if size(t,2)==1
    cond2a=(pfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 0],[0 0 kappafunc(s(1),kappaccp)],s(1),k1,k2,sig1,sig2)*x0'+qfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 0],[0 0 kappafunc(s(1),kappaccp)],s(1),k1,k2,sig1,sig2,theta1,theta2))*exp(bfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 0],s(1),k1,k2,sig1,sig2)*x0'+cfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 0],s(1),k1,k2,sig1,sig2,theta1,theta2))*(1-w20*etafunc(s(1),expectedpsiccp));
    cond2b=(pfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 -1],[0 0 kappafunc(s(1),kappaccp)],s(1),k1,k2,sig1,sig2)*x0'+qfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 -1],[0 0 kappafunc(s(1),kappaccp)],s(1),k1,k2,sig1,sig2,theta1,theta2))*exp(bfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 -1],s(1),k1,k2,sig1,sig2)*x0'+cfunc(0,v,[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],[0 0 -1],s(1),k1,k2,sig1,sig2,theta1,theta2))*w21*etafunc(s(1),expectedpsiccp);
else
    nn=size(t,2);%number of intervals
    ba=bfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 0],s(nn),k1,k2,sig1,sig2)*h(s(nn-1),s(nn),rtimes,l1times,l2times);
    pa=pfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 0],[0 0 kappafunc(s(nn),kappaccp)],s(nn),k1,k2,sig1,sig2)*h(s(nn-1),s(nn),rtimes,l1times,l2times);
    cas=zeros(1,nn);
    qas=zeros(1,nn);
    cas(nn)=cfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 0],s(nn),k1,k2,sig1,sig2,theta1,theta2);
    qas(nn)=qfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 0],[0 0 kappafunc(s(nn),kappaccp)],s(nn),k1,k2,sig1,sig2,theta1,theta2);

    bb=bfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 -1],s(nn),k1,k2,sig1,sig2)*h(s(nn-1),s(nn),rtimes,l1times,l2times);
    pb=pfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 -1],[0 0 kappafunc(s(nn),kappaccp)],s(nn),k1,k2,sig1,sig2)*h(s(nn-1),s(nn),rtimes,l1times,l2times);
    cbs=zeros(1,nn);
    qbs=zeros(1,nn);
    cbs(nn)=cfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 -1],s(nn),k1,k2,sig1,sig2,theta1,theta2);
    qbs(nn)=qfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) kappafunc(s(nn),kappaccp)],[0 0 -1],[0 0 kappafunc(s(nn),kappaccp)],s(nn),k1,k2,sig1,sig2,theta1,theta2);
    for n=(nn-1):-1:2
        cas(n)=cfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],ba,s(n),k1,k2,sig1,sig2,theta1,theta2);
        qas(n)=qfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],ba,pa,s(n),k1,k2,sig1,sig2,theta1,theta2);
        pa=pfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],ba,pa,s(n),k1,k2,sig1,sig2)*h(s(n-1),s(n),rtimes,l1times,l2times);
        ba=bfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],ba,s(n),k1,k2,sig1,sig2)*h(s(n-1),s(n),rtimes,l1times,l2times);
        
        cbs(n)=cfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],bb,s(n),k1,k2,sig1,sig2,theta1,theta2);
        qbs(n)=qfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],bb,pb,s(n),k1,k2,sig1,sig2,theta1,theta2);
        pb=pfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],bb,pb,s(n),k1,k2,sig1,sig2)*h(s(n-1),s(n),rtimes,l1times,l2times);
        bb=bfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) kappafunc(s(n),kappaccp)],bb,s(n),k1,k2,sig1,sig2)*h(s(n-1),s(n),rtimes,l1times,l2times);
    end
    cas(1)=cfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],ba,s(1),k1,k2,sig1,sig2,theta1,theta2);
    qas(1)=qfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],ba,pa,s(1),k1,k2,sig1,sig2,theta1,theta2);
    pa=pfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],ba,pa,s(1),k1,k2,sig1,sig2);
    ba=bfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],ba,s(1),k1,k2,sig1,sig2);
    
    cbs(1)=cfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],bb,s(1),k1,k2,sig1,sig2,theta1,theta2);
    qbs(1)=qfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],bb,pb,s(1),k1,k2,sig1,sig2,theta1,theta2);
    pb=pfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],bb,pb,s(1),k1,k2,sig1,sig2);
    bb=bfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) kappafunc(s(1),kappaccp)],bb,s(1),k1,k2,sig1,sig2);
    
    cond2a=(sum(qas)+pa*x0')*exp(ba*x0'+sum(cas))*(1-w20*etafunc(s(nn),expectedpsiccp));
    cond2b=(sum(qbs)+pb*x0')*exp(bb*x0'+sum(cbs))*w21*etafunc(s(nn),expectedpsiccp);

 
end