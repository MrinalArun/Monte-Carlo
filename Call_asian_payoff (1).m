function h=Call_asian_payoff(Z,r,sigma,T,m,S0,K)
% payoff of and Asian Call option
A = (r-0.5*sigma^2)*T/m*transpose(1:m)+sigma*sqrt(T/m)*cumsum(Z);
g = mean(exp(A))*S0;
h = max(0,g-K);
end