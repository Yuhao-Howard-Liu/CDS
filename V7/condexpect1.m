function [cond1] = condexpect1(v,t,s,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,rtimes,l1times,l2times)
%for (6.1) in v38
%v scalar -- the discretized time points
%t vector -- the times for the MC relization
%s vector -- the states for the MC relization
t=t(t<=v);
s=s(t<=v);

if size(t,2)==1
    cond1=exp(bfunc(0,v,[1 kappafunc(s(1),kapparf) 0],[0 0 0],s(1),k1,k2,sig1,sig2)*x0'+cfunc(0,v,[1 kappafunc(s(1),kapparf) 0],[0 0 0],s(1),k1,k2,sig1,sig2,theta1,theta2));
else
    nn=size(t,2);%number of intervals
    b=bfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) 0],[0 0 0],s(nn),k1,k2,sig1,sig2)*h(s(nn-1),s(nn),rtimes,l1times,l2times);
    cs=zeros(1,nn);
    cs(nn)=cfunc(t(nn),v,[1 kappafunc(s(nn),kapparf) 0],[0 0 0],s(nn),k1,k2,sig1,sig2,theta1,theta2);

    for n=(nn-1):-1:2
        cs(n)=cfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) 0],b,s(n),k1,k2,sig1,sig2,theta1,theta2);
        b=bfunc(t(n),t(n+1),[1 kappafunc(s(n),kapparf) 0],b,s(n),k1,k2,sig1,sig2)*h(s(n-1),s(n),rtimes,l1times,l2times);
    end
    cs(1)=cfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) 0],b,s(1),k1,k2,sig1,sig2,theta1,theta2);
    b=bfunc(t(1),t(2),[1 kappafunc(s(1),kapparf) 0],b,s(1),k1,k2,sig1,sig2);
    cond1=exp(b*x0'+sum(cs));

 
end