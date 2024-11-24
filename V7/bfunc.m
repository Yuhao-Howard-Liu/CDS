function [bb] = bfunc(t,T,alp,u,env,k1,k2,sig1,sig2)
    % function \boldsymbol{b}
    u = u(:)';
    k = kfunc(env,k1,k2);
    sig = sigmafunc(env,sig1,sig2);
    
    h = sqrt(2 * sig.^2 .* alp + k.^2);
    a = (h - (-sig.^2 .* u + k)) ./ (h + (-sig.^2 .* u + k));
    
    exp_term = exp(-h .* (T - t));
    
    bb = k ./ sig.^2 - (h .* (1 - exp_term .* a)) ./ (sig.^2 .* (1 + exp_term .* a));
end
