function [R, LS, LI, At] = MVA_exact(lambda, pij, Tij, m)
% build MOD network
% vector of size n for the SS nodes
% matrix of size n x n for the IS nodes
% lambda = service rates for SS nodes
% 1/Tij = service rates for IS nodes
% each node requires: 
% Response time (or waiting time) 
%   - only need response time for SS nodes. For IS nodes, it's just Tij
% Mean queue length
n = length(lambda);

viss = null((eye(n) - pij'));
viss = viss./viss(1);
viis2 = pij.*repmat(viss, 1, n);

% R = response times for SS nodes
% Tij = response times for IS nodes
% LS = queue lengths for SS nodes
% LI = queue lengths for IS nodes
% At = availability for SS nodes = (LS/R)/mu
R = zeros(n, m);
LS = zeros(n, m);
LI = zeros(n,n);
LITotal = zeros(m,1);
At = zeros(n, m);
for i = 1:m
    
    if i == 1
        R(:,i) = 1./lambda;
    else
        R(:,i) = 1./lambda + LS(:,i-1)./lambda;
    end
    Lden = sum(viss.*R(:,i)) + sum(sum(viis2.*Tij));
    LS(:,i) = i.*(viss.*R(:,i))/Lden;
    LI = i/Lden.*(viis2.*Tij);
    LITotal(i) = sum(sum(LI));
    i;
    At(:,i) = (LS(:,i)./R(:,i))./lambda;
    
end

end