function [ss] = sigmafunc(env,sig1,sig2)
% second diagonal elements of Sigma
st = env == 1;
ss = sig1.*st + sig2.*(1-st);
% if state==1
%     ss=sig1;
% elseif state ==2
%     ss=sig2;
% end
end