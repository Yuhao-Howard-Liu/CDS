function [tt] = thetafunc(env,theta1,theta2)
st = env == 1;
tt = theta1.*st + theta2.*(1-st);
% if state==1
%     tt=theta1;
% elseif state ==2
%     tt=theta2;
% end
end