function [EstPrices] = GBM_SDE_KnockoutPut()
% This code implements an Euler SDE scheme for the pricing of a European knockout put option
% in the Black-Scholes framework. Therefore have dX_t = r X_t dt + sig X_t dW_t
% Of course we know the solution of the SDE so we can compare the true knockout put price from 
% Black-Scholes with the estimated price provided by (i) the usual Euler scheme and (ii) the Euler 
% scheme where we also simulate the minimum on the price process on each interval using equation (6.50)
% in Glasserman (With a "-" instead of "+" sign in the appropriate position.)

% Input Parameters:  X0 = initial stock price, T = maturity, K = strike, sig = vol, r = riskfree rate = drift
X0 = 100;       K = 90;    T = .5;   q = 0;   r = .01;   sig = .5;  H = 80;  % Must have barrier H < strike K
%randn('state',sum(100*clock));
%rng('default') ; % Sets random number generator seed to default value
TruePutKnockOut = BS_KnockoutPut(X0,sig,r,q,T,K,H)
NumSteps = [10 50 100 250 500 1000];
dt = T./NumSteps;
EstPrices = zeros(2,length(NumSteps));
NumSamples = 4000000;  % Number of samples for each dt. Should take around 2m

for i=1:length(NumSteps)
    [EulerEst,~,EulerMinEst,~] = EulerScheme(X0,sig,r,q,dt(i),NumSamples,NumSteps(i),K,T,H);
    EstPrices(1,i) = EulerEst;
    EstPrices(2,i) = EulerMinEst;
    i
end
Abs_Error = abs(TruePutKnockOut - EstPrices);

% Now plot the absolute error for the schemes as a function of NumSteps
for i=1:2
    loglog(NumSteps,Abs_Error(i,:),'linewidth',1);
    hold on
end
grid on

set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Steps','fontsize',15)
ylabel('Absolute Error','rot',90,'fontsize',15)
legend('Standard Euler','Euler-Min');
end


%%**************************************************************************************
function [EulerEst,EulerStd,EulerMinEst,EulerMinStd] = EulerScheme(X,sig,r,q,dt,n,NumSteps,K,T,H)

Con1 = sig*sqrt(dt);
Xold = X;
MinPrice = X;
BetterMinPrice = X;
for j=1:NumSteps  % NumPeriods corresponds to l
    Xnew = Xold.*(1 + (r-q)*dt + Con1*randn(n,1));    
    Xold = Xnew;
    MinPrice = min(MinPrice, Xnew);
    LocalMin = (Xnew + Xold - sqrt((Xnew - Xold).^2 - 2*dt*((sig*Xold).^2).*log(rand(n,1)))) / 2;
    BetterMinPrice = min(BetterMinPrice, LocalMin);
end

DiscountValue = exp(-r*T)*max(0,K - Xnew) .* [MinPrice > H];
EulerEst = mean(DiscountValue);
EulerStd = std(DiscountValue)/sqrt(n);
%EulerCI = [EulerEst - 1.96*EulerStd   EulerEst + 1.96*EulerStd];

DiscountValue = exp(-r*T)*max(0,K - Xnew) .* [BetterMinPrice > H];  % Now use superior min estimate
EulerMinEst = mean(DiscountValue);
EulerMinStd = std(DiscountValue)/sqrt(n);
end

%%**************************************************************************************
function [PutKnockOut] = BS_KnockoutPut(S,sigma,r,q,T,K,H)
% This function prices a knockout (continuously sampled) put option with strike = K, barrier = H and 
% time-to-maturity = T. Current spot price = S. The formula is taken as p_do = p - p_di where p is the
% value of a vanilla put, and p_di is the value of the corresponding down-and-in put option. The formula for
% p_di is taken from Hull. H < K is assumed otherwise price = 0.

d1 = (log(S/K) + (r - q + sigma.^2/2)*T) / (sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);

lambda = (r-q+sigma.^2/2) / (sigma.^2);
y = log((H^2)/(S*K)) / (sigma*sqrt(T))  +  lambda*sigma*sqrt(T);
x1 = log(S/H) / (sigma*sqrt(T))  +  lambda*sigma*sqrt(T);
y1 = log(H/S) / (sigma*sqrt(T))  +  lambda.*sigma*sqrt(T);

PutPrice = K*exp(-r*T)*normcdf(-d2) - S*exp(-q*T)*normcdf(-d1);

PutKnockIn = -S*normcdf(-x1)*exp(-q*T) + K*exp(-r*T)*normcdf(-x1+sigma*sqrt(T)) + ...
    S*exp(-q*T)*((H/S)^(2*lambda))*(normcdf(y) - normcdf(y1)) - ...
    K*exp(-r*T)*((H/S)^(2*lambda-2))*(normcdf(y-sigma*sqrt(T)) - normcdf(y1-sigma*sqrt(T)));

PutKnockOut = PutPrice - PutKnockIn;

end