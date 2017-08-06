function [] = GBM_SDE_Schemes()
% This code implements various SDE schemes for the pricing of a European call option in
% in the Black-Scholes framwork. Therefore have dX_t = r X_t dt + sig X_t dW_t
% Of course we know the solution of the SDE so we can compare the true call option price from 
% Black-Scholes with the estimated price provided by (i) the usual Euler scheme and (ii) the Euler 
% scheme with Richardson Extrapolation 

% Input Parameters:  X0 = initial stock price, T = maturity, K = strike, sig = vol, r = riskfree rate = drift
X0 = 100;       K = 100;    T = .5;      r = .01;    sig = .4;
%randn('state',sum(100*clock));
%rng('default') ; % Sets random number generator seed to default value
TrueCallPrice =  BlackScholes(X0,sig,r,T,K)
NumSteps = [10 50 100 250 500];
dt = T./NumSteps;
EstPrices = zeros(2,length(NumSteps));
NumSamples = 16000000;  % Number of samples for each dt. Should take around 2m

for i=1:length(NumSteps)
    [EstPrice, StdPrice, CI] = EulerEst(X0,sig,r,dt(i),NumSamples,NumSteps(i),K,T);
    EstPrices(1,i) = EstPrice;
    
    [EstRichPrice, StdRichPrice, CI_Rich] = EulerRichardson(X0,sig,r,dt(i),NumSamples,NumSteps(i),K,T);
    EstPrices(2,i) = EstRichPrice;
    i
end
Abs_Error = abs(TrueCallPrice - EstPrices);

% Now plot the absolute error for the schemes as a function of NumSteps
for i=1:2
    loglog(NumSteps,Abs_Error(i,:),'linewidth',1);
    hold on
end
grid on

set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Steps','fontsize',15)
ylabel('Absolute Error','rot',90,'fontsize',15)
legend('Standard Euler','Euler-Richardson');
end

%***************************************************************************
function [Rich_Est,Rich_std,Rich_CI] = EulerRichardson(X,sig,r,dt,n,NumSteps,K,T)

Z = X;  % for the 2*dt calculations
Con1 = sig*sqrt(dt);

for j=1:NumSteps  
    New_BM_Inc = Con1*randn(n,1);
    X = X.*(1 + r*dt + New_BM_Inc);
    
    if (mod(j,2) == 0)
        Z = Z.*(1 + 2*r*dt + Old_BM_Inc + New_BM_Inc);
    end
    Old_BM_Inc = New_BM_Inc;
end
Rich_Samples = exp(-r*T)*(2*max(0,X-K) - max(0, Z-K));
Rich_Est = mean(Rich_Samples);    %exp(-r*T)*(2*mean(max(0,X-K)) - mean(max(0, Z-K)));
Rich_std = std(Rich_Samples)/sqrt(n);
Rich_CI = [Rich_Est - 1.96*Rich_std   Rich_Est + 1.96*Rich_std];
end

%%**************************************************************************************
function [EulerEst,EulerStd,EulerCI] = EulerEst(X,sig,r,dt,n,NumSteps,K,T)

Con1 = sig*sqrt(dt);
for j=1:NumSteps  % NumPeriods corresponds to l
    X = X.*(1 + r*dt + Con1*randn(n,1));
end

DiscountValue = exp(-r*T)*max(0,X-K);
EulerEst = mean(DiscountValue);
EulerStd = std(DiscountValue)/sqrt(n);
EulerCI = [EulerEst - 1.96*EulerStd   EulerEst + 1.96*EulerStd];
end

%%**************************************************************************************
function [CallPrice] = BlackScholes(S,sigma,r,T,K)

% S and T are vectors
q = 0; % Div yield set to 0 throughout these experiments
d1 = (log(S./K) + (r-q + (sigma.^2)/2).*T) ./ (sigma.*sqrt(T));
d2 = d1 - sigma.*sqrt(T);
CallPrice = (S.*exp(-q.*T).*normcdf(d1) - K.*exp(-r.*T).*normcdf(d2));

end