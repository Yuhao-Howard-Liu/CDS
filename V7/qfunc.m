function [ccc] = qfunc(t,T,alp,u,v,env,k1,k2,sig1,sig2,theta1,theta2)
    % function q
    v = v(:)';
    u = u(:)';
    k = kfunc(env,k1,k2);
    sig = sigmafunc(env,sig1,sig2);
    theta = thetafunc(env,theta1,theta2);
    
    h = sqrt(2 * sig.^2 .* alp + k.^2);
    a = (h - (-sig.^2 .* u + k)) ./ (h + (-sig.^2 .* u + k));
    
    term1 = 4 * (k .* theta .* v ./ sig.^2) .* h ./ (2 * alp - (sig .* u).^2 + 2 * u .* k) ./ (1 + exp(-h * (T - t)) .* a);
    term2 = 2 * (k .* theta .* v ./ sig.^2) .* (h + (-sig.^2 .* u + k)) ./ (2 * alp - (sig .* u).^2 + 2 * u .* k);
    
    ccc = sum(term1 - term2);
end
