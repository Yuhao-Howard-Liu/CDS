function [bbb] = pfunc(t,T,alp,u,v,env,k1,k2,sig1,sig2)
    % function p
    v = v(:)';
    u = u(:)';
    k = kfunc(env,k1,k2);
    sig = sigmafunc(env,sig1,sig2);
    
    h = sqrt(2 * sig.^2 .* alp + k.^2);
    a = (h - (-sig.^2 .* u + k)) ./ (h + (-sig.^2 .* u + k));
    exp_term = exp(-h .* (T - t));
    
    bbb = v .* exp_term .* ((1 + a) ./ (1 + a .* exp_term)).^2;
end