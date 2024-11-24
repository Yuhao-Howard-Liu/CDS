function [gg] = etafunc(env,kk)
% expected value of expected value of e^-Psi in different states
st = env == 1;
gg = kk(1).*st + kk(2).*(1-st);
end