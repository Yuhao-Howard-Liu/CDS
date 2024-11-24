function [t,s] = ctmcgenerator(tmax,initial_s,Q)
t=0;
s=initial_s;%initial state
tend=t(end);
while tend<tmax
    i=s(end);
    q=-Q(i,i);
    if q==0
        break
    else
       tend =-log(rand)/(-Q(i,i))+t(end); %exponential holding time
       if tend<tmax
       t=[t tend];
       p=Q(i,:);
       p(i)=0;
       p=p./sum(p);
       s=[s samplefromp(p,1)];
       end
    end
end
end

