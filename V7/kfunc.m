function [gg] = kfunc(env,k1,k2)
% diagonal elements of K
st = env == 1;
gg = k1.*st + k2.*(1-st);
end

