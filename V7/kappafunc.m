function [gg] = kappafunc(env,kk)
st = env == 1;
gg = kk(1).*st + kk(2).*(1-st);
% kk=kk(:);
% if state==1
%     gg=kk(1);
% elseif state ==2
%     gg=kk(2);
end