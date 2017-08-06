% This script implements the combination of pathwise and LR estimators for estimating the delta 
% of a digital European call option with strike K in the Black-Scholes framework. The estimator
% takes the form (1/2*)eps * 1_{|S_T -K| < eps}*S_T / S_0  + h_{eps}(S_T) zeta(S_T)/(S_0*sig*sqrt(T)) 
% with  h_{eps}(x) := 1_{x > K} - min{1, max{0, x-K+eps} / (2*eps) and 
% zeta(x) := (log(x/S_0) -(r-sig^2/2)*T) / (sig*sqrt(T)). See page 417 of Glasserman

% Input Parameters:  X0 = initial stock price, T = maturity, K = strike, sig = vol, r = riskfree rate = drift
S0 = 100;   K = 100;   r = .05;    sig = .3;    T = 0.25;

rng('default') ; % Sets random number generator seed to default value. Can also randomize it via "randn('state',sum(100*clock));"
eps = 10:70;  % Range of epsilon values
NumSamples = 100000;
Con1 = (r-(sig^2)/2)*T; 

ST = S0*exp(Con1 + sig*sqrt(T)*randn(NumSamples,1));
Zeta = (log(ST/S0) - (r-sig^2/2)*T) / (sig*sqrt(T));
Pathwise_Term = abs(ST - K) .* ST/(2*S0);
LR_Term = Zeta / (S0*sig*sqrt(T));
h_Term1 = [ST > K];
Delta_Estimators = zeros(1,length(eps));   Delta_Estimator_Vars = zeros(1,length(eps)); 

for i=1:length(eps)
    MixedEstimator = exp(-r*T)* ([Pathwise_Term < eps(i)]/ eps(i)) ...
        + LR_Term .* (h_Term1 - min(1, max(0, ST-K+eps(i)) / (2*eps(i))));
    Delta_Estimators(i) = mean(MixedEstimator);
    Delta_Estimator_Vars(i) = var(MixedEstimator);
end

% Now produce Fig 7.4 of Glasserman 
plot(eps,Delta_Estimator_Vars,'linewidth',1);
hold on
grid on
set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Epsilon','fontsize',15)
ylabel('Variance','rot',90,'fontsize',15)
