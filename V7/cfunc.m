function [cc] = cfunc(t,T,alp,u,env,k1,k2,sig1,sig2,theta1,theta2)
    % function c
    u = u(:)';
    k = kfunc(env,k1,k2);
    sig = sigmafunc(env,sig1,sig2);
    theta = thetafunc(env,theta1,theta2);
    
    h = sqrt(2 * sig.^2 .* alp + k.^2);
    
    a = (h - (-sig.^2 .* u + k)) ./ (h + (-sig.^2 .* u + k));
    
    exp_term = exp(-h .* (T - t));
    
    cc = - sum(k .* theta ./ sig.^2 .* (h - k) .* (T - t) + 2 * k .* theta ./ sig.^2 .* (log(1 + exp_term .* a) - log(1 + a)));
end
