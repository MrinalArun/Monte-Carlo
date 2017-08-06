clc;
r = 0.05;
sigma = 0.25;
T = 1;
m = 6;
S0 = 100;
K = 90;
N = 1000;   % samples for each stratum
N_s = 10;   % number of intervals for W

% calculate weight Mu
A = exp((r-sigma^2*0.5)*T/m); 
B = sigma*sqrt(T/m);
Temp = zeros(m,1); 
Mu = zeros(m,1);
u = Temp;
for i=1:m
    Temp(i)=A^(m+1-i);
end
Temp = cumsum(Temp);
for i=1:m
    Mu(i) = B*S0*Temp(m+1-i)/m;
end
Mu = 1/norm(Mu)*Mu;

U = 0:1/N_s:1;
H = zeros(N,N_s);
for i=1:N_s
    a = U(i); %left bound  of strata
    W = norminv(rand(N,1)*(1/N_s)+a); %generate W
    for j=1:N
        V = randn(m,1);
        Z = W(j)*Mu + V - Mu*Mu'*V; %generate X
        H(j,i) = Call_asian_payoff(Z,r,sigma,T,m,S0,K);
    end
end
price = mean(mean(H))*exp(-r*T);
price_stderr =sqrt(mean(std(H).^2)/N)*exp(-r*T);
alpha = norminv(0.99);
disp('simulation result for stratified sampling: ');
disp(['estimation of price is ',num2str(price)]);
disp(['std error of price is ',num2str(price_stderr)]);
disp(['99% Confidence interval is [ ',num2str(price-alpha*price_stderr),' , ',num2str(price+alpha*price_stderr),' ] .']);
      