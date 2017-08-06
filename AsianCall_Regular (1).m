clc;
% Q5(e) regular sampling for asian call option pricing
r = 0.05;
sigma = 0.25;
T = 1;
m = 6;
S0 = 100;
K = 90;
N = 19000;   % sample size

%% Regular Sampling
U = 0:1/N_s:1;
H = zeros(N,1);
for j=1:N
    Z = randn(m,1);
    H(j,1) = Call_asian_payoff(Z,r,sigma,T,m,S0,K);
end
price = mean(H)*exp(-r*T);
price_stderr = std(H)/sqrt(N)*exp(-r*T);
alpha = norminv(0.99);
disp('simulation result for regular sampling: ');
disp(['estimation of price is ',num2str(price)]);
disp(['standard error of price is ',num2str(price_stderr)]);
disp(['99% Confidence interval is [ ',num2str(price-alpha*price_stderr),' , ',num2str(price+alpha*price_stderr),' ] .']);

        
        
        
        
        
        