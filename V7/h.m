function [vv] = h(s1,s2,rtimes,l1times,l2times)
%The diagonal vector of h for different s1 and s2, rtimes is the
%multiples between good and bad states for risk free rate , ltimes is the
%multiples between good and bad states for  lamabda
if s1==1 
    if s2==1
        vv=diag([1 1 1]);
    elseif s2==2
        vv=diag([rtimes l1times l2times]);
    end
elseif s1==2
    if s2==1
        vv=diag([1/rtimes 1/l1times 1/l2times]);
    elseif s2==2
        vv=diag([1 1 1]);
    end
end