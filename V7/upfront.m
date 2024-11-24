function [premium,CCPloss,protection] = upfront(M,N,Tmax,initialstate,Q,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsirf,expectedpsiccp,w10,w11,w20,w21)


% Preallocate arrays for efficiency
NN=N+1;
PS1 = zeros( NN, M);
PS2 = zeros( NN, M);
PS3 = zeros( NN, M);

% Start parallel processing
parfor m = 1:M
    % Generate continuous-time Markov chain sample paths
    [t, s] = ctmcgenerator(Tmax, initialstate, Q);
    % Define the time intervals outside the loop to avoid recomputation
    timeIntervals = [0 (1:N) * Tmax / N];
    for n = 1:NN
        [cond1] = condexpect1(timeIntervals(n),t,s,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,rtimes,l1times,l2times)
        [cond2a, cond2b] = condexpect2(timeIntervals(n),t,s,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,kappaccp,rtimes,l1times,l2times,expectedpsiccp,w20,w21)
        [cond5a, cond5b] = condexpect5(timeIntervals(n),t,s,x0,k1,k2,sig1,sig2,theta1,theta2,kapparf,rtimes,l1times,l2times,expectedpsirf,w10,w11)
        % Store results in preallocated matrices
        PS1(n, m) = cond1; 
        PS2(n, m) = cond2a-cond2b;
        PS3(n, m) = cond5a-cond5b;
    end
end



PS1=mean(PS1,2);
PS2=mean(PS2,2);
PS3=mean(PS3,2);
timeIntervals = [0 (1:N) * Tmax / N];
premium=trapz(timeIntervals,PS1);
CCPloss=trapz(timeIntervals,PS2);
protection=trapz(timeIntervals,PS3);








end