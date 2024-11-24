sppaths=[1 2 1];
sppatht=[0 0.7 1.5];
x0=[0.05 0.02 0.001];
T=2;
N=41;
dt=T/(N-1);
sig1=[0.04 0.15 0.005];
sig2=[0.06 0.18 0.02];
theta1_p=[0.05 0.02 0.001];
theta2_p=[0.025 0.04 0.005];
k1_p=[0.1 1 1];
k2_p=[0.12 2 2];
% timepoints=0:dt:2;
X = zeros(N, 3);
X(1,:) = x0;
    
for i = 1:13
    dW = sqrt(dt) * randn(1,3);  % Brownian increment
    % Ensure the square root is real
    X(i+1,1) = X(i,1) + k1_p(1)*(theta1_p(1) - max(X(i,1), 0))*dt + sig1(1)*sqrt(max(X(i,1), 0))*dW(1);
    X(i+1,2) = X(i,2) + k1_p(2)*(theta1_p(2) - max(X(i,2), 0))*dt + sig1(2)*sqrt(max(X(i,2), 0))*dW(2);
    X(i+1,3) = X(i,3) + k1_p(3)*(theta1_p(3) - max(X(i,3), 0))*dt + sig1(3)*sqrt(max(X(i,3), 0))*dW(3);
end

dW = sqrt(dt) * randn(1,3);
X(15,1) = X(14,1)/2 + k1_p(1)*(theta1_p(1) - max(X(14,1), 0))*dt + sig1(1)*sqrt(max(X(14,1), 0))*dW(1);
X(15,2) = X(14,2)*2 + k1_p(2)*(theta1_p(2) - max(X(14,2), 0))*dt + sig1(2)*sqrt(max(X(14,2), 0))*dW(2);
X(15,3) = X(14,3)*5 + k1_p(3)*(theta1_p(3) - max(X(14,3), 0))*dt + sig1(3)*sqrt(max(X(14,3), 0))*dW(3);

for i = 15:29
    dW = sqrt(dt) * randn(1,3);  % Brownian increment
    % Ensure the square root is real
    X(i+1,1) = X(i,1) + k2_p(1)*(theta2_p(1) - max(X(i,1), 0))*dt + sig2(1)*sqrt(max(X(i,1), 0))*dW(1);
    X(i+1,2) = X(i,2) + k2_p(2)*(theta2_p(2) - max(X(i,2), 0))*dt + sig2(2)*sqrt(max(X(i,2), 0))*dW(2);
    X(i+1,3) = X(i,3) + k2_p(3)*(theta2_p(3) - max(X(i,3), 0))*dt + sig2(3)*sqrt(max(X(i,3), 0))*dW(3);
end

dW = sqrt(dt) * randn(1,3);  % Brownian increment
% Ensure the square root is real
X(31,1) = X(30,1)*2 + k2_p(1)*(theta2_p(1) - max(X(30,1), 0))*dt + sig2(1)*sqrt(max(X(30,1), 0))*dW(1);
X(31,2) = X(30,2)/2 + k2_p(2)*(theta2_p(2) - max(X(30,2), 0))*dt + sig2(2)*sqrt(max(X(30,2), 0))*dW(2);
X(31,3) = X(30,3)/5 + k2_p(3)*(theta2_p(3) - max(X(30,3), 0))*dt + sig2(3)*sqrt(max(X(30,3), 0))*dW(3);

for i = 31:40
    dW = sqrt(dt) * randn(1,3);  % Brownian increment
    % Ensure the square root is real
    X(i+1,1) = X(i,1) + k1_p(1)*(theta1_p(1) - max(X(i,1), 0))*dt + sig1(1)*sqrt(max(X(i,1), 0))*dW(1);
    X(i+1,2) = X(i,2) + k1_p(2)*(theta1_p(2) - max(X(i,2), 0))*dt + sig1(2)*sqrt(max(X(i,2), 0))*dW(2);
    X(i+1,3) = X(i,3) + k1_p(3)*(theta1_p(3) - max(X(i,3), 0))*dt + sig1(3)*sqrt(max(X(i,3), 0))*dW(3);
end
plot(0:dt:T, X);

save xzero.mat X