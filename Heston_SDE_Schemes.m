function [] = Heston_SDE_Schemes()
% This code implements various SDE schemes for the pricing of a European call option in
% in the Heston model with  dS_t / S_t = (r - q) dt + sigma_t dW_t   and 
% d sigma_t^2 = kappa * (eta - sigma_t^2) * dt + theta * sigma_t dW^~_t  
% where the instantaneous correlation between the 2 Q-Brownian motions = rho.
% We can compare the true call option price from Fourier inversion of characteristic function of the log stock
% price with the estimated price provided by (i) the usual Euler scheme and (ii) the second order 
% scheme provided in e.g. 6.2.2 of Glasserman

% randn('state',sum(100*clock));
%rng('default') ; % Sets random number generator seed to default value

% Input Parameters:  
global S0 K r T q sigma_0 theta rho kappa eta  
%************************************************************************
% Params on next line are for Glasserman e.g. 6.2.2
S0 = 100;  K = 100;   r = .05;     T = 1;   q = 0;   
sigma_0 = sqrt(.04);    theta = .3;       rho = -.5;     kappa = 1.2;      eta = .04;  
% Parameters on next 2 lines are for Andersen's 2007 example
S0 = 100; K = 100;  T = 10; r = 0; q = 0; 
sigma_0 = sqrt(.04);    theta = .1;       rho = -.9;     kappa = .5;      eta = .04; 
%**************************************************************************
alp = .75;  % This is the "integration" parameter 

F = @(v) (exp(-r*T - 1i*v*log(K))) .* HestonCharFunction(v-(alp+1)*1i)  ./ (alp^2 + alp - v.^2 + 1i*(2*alp+1)*v);
TrueCallPrice = (exp(-alp*log(K))/pi) * real(quad(F,0,1000)) % Need real part 

NumSteps = [4 8 16 30 100 500 1000];    % Elements in NumSteps should be even so that Richardson extrapolation works 
dt = T./NumSteps;
EstPrices = zeros(2,length(NumSteps));
NumSamples = 1000000;  % Number of samples for each dt. Should take around 2m

for i=1:length(NumSteps)
    EstPrices(1,i) = Euler_Heston(NumSamples,dt(i),NumSteps(i));
    EstPrices(2,i) = Euler_Richardson(NumSamples,dt(i),NumSteps(i));
    EstPrices(3,i) = SecondOrder(NumSamples,dt(i),NumSteps(i));
    i
end
Abs_Error = abs(TrueCallPrice - EstPrices);

% Now plot the absolute error for the schemes as a function of NumSteps
for i=1:3
    loglog(NumSteps,Abs_Error(i,:),'linewidth',1);
    hold on
end
grid on

set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Steps','fontsize',15)
ylabel('Absolute Error','rot',90,'fontsize',15)
legend('Euler','Euler-Richardson','Second Order');
end

%%**************************************************************************************
function [SecondOrderEst] = SecondOrder(NumSamples,dt,NumSteps)
% This is second-order scheme from E.G. 6.2.2 of Glasserman's Monte Carlo text
global S0 K r T q sigma_0 theta rho kappa eta  

X = S0;     Var = sigma_0^2;     C1 = sqrt(1 - rho^2);
sig1 = theta*rho;   sig2 = theta*C1;
C2 = ((r-q)*dt)^2/2;    C3 = (r-q +(sig1-kappa)/4);   C4 = kappa*eta/4 - theta^2/16; 

for i=1:NumSteps
    Z = randn(NumSamples,2);
    Tm1 = sqrt(Var * dt);
    Tm2 = sqrt(dt./Var);
    
    Term1 = 1 + (r - q) * dt + Tm1 .* Z(:,1) + C2;
    Term2 = (C3*Tm1 + C4*Tm2).*Z(:,1)*dt;
    UniformPiece = 1 - 2*[rand(NumSamples,1) <= .5];
    Term3 = (Var+sig1/2).*(Z(:,1).^2 - 1)*dt/2  + (sig2/4)*dt*(Z(:,1).*Z(:,2) + UniformPiece);
    X = X.*(Term1 + Term2 + Term3);
    
    Var = Var + kappa*(eta - Var) * dt + theta * Tm1 .* (rho* Z(:,1) + C1 * Z(:,2)) - ...
        (kappa^2)*(eta-Var)*(dt^2)/2 + (C4*Tm2 - (3*kappa/2)*Tm1) * theta .* (rho* Z(:,1) + C1 * Z(:,2))*dt ...
        + (sig1^2/4)*dt*(Z(:,1).^2-1) + (sig2^2/4)*dt*(Z(:,2).^2-1) + ((sig1*sig2)/2)*dt*Z(:,1).*Z(:,2);    
    Var = abs(Var);
end

DiscountValue = exp(-r*T)*max(0,X-K);
SecondOrderEst = mean(DiscountValue);

end

%***************************************************************************
function [Richardson_Est] = Euler_Richardson(NumSamples,dt,NumSteps)

global S0 K r T q sigma_0 theta rho kappa eta  
X = S0;     Var = sigma_0^2;        Con1 = sqrt(1 - rho^2);
X_2dt = X;    Var_2dt = Var;
for i=1:NumSteps
    New_BM_Inc = randn(NumSamples,2);
    
    Term1 = sqrt(Var * dt);
    X = X.*(1 + (r - q) * dt + Term1 .* New_BM_Inc(:,1));
    Var = Var + kappa*(eta - Var) * dt + theta * Term1 .* (rho* New_BM_Inc(:,1) + Con1 * New_BM_Inc(:,2));    
    Var = abs(Var);
    
    if (mod(i,2) == 0)
        Term1_2dt = sqrt(Var_2dt * dt);
        X_2dt = X_2dt.*(1 + 2*(r - q) * dt + Term1_2dt .* (Old_BM_Inc(:,1) + New_BM_Inc(:,1)));
        
        Var_2dt = Var_2dt + kappa*(eta - Var_2dt) * 2*dt + theta * Term1_2dt .* ...
            (rho* (Old_BM_Inc(:,1) + New_BM_Inc(:,1)) + Con1 * (Old_BM_Inc(:,2) + New_BM_Inc(:,2))); 
        Var_2dt = abs(Var_2dt);
    end
    Old_BM_Inc =  New_BM_Inc;
end

Richardson_Est = exp(-r*T)*(2*mean(max(0,X-K)) - mean(max(0, X_2dt-K)));
end

%%**************************************************************************************
function [EulerEst] = Euler_Heston(NumSamples,dt,NumSteps)

global S0 K r T q sigma_0 theta rho kappa eta  
X = S0;     Var = sigma_0^2;        Con1 = sqrt(1 - rho^2);

for i=1:NumSteps
    Z = randn(NumSamples,2);
    Term1 = sqrt(Var * dt);
    X = X.*(1 + (r - q) * dt + Term1 .* Z(:,1));
    Var = Var + kappa*(eta - Var) * dt + theta * Term1 .* (rho* Z(:,1) + Con1 * Z(:,2));    
    Var = abs(Var);
end

EulerEst = mean(exp(-r*T)*max(0,X-K));
end

%%**************************************************************************************
function [CharFun] = HestonCharFunction(u)
     
global S0 r T q sigma_0 kappa eta rho theta
  
d = sqrt((rho*theta*u*1i - kappa).^2 - (theta^2)*(-1i*u - u.^2));
g = (kappa - rho*theta*u*1i-d) ./ (kappa - rho*theta*u*1i +d);
% CharFun = E[exp(iu log(S_t) | S_0, sigma_0)] for vector u.
CharFun = exp(1i*u*(log(S0) + (r-q)*T) + ...
    eta*kappa*(theta^(-2))*((kappa- rho*theta*u*1i -d)*T  - 2*log((1-g.*exp(-d*T))./(1-g))) + ...
    sigma_0^2 * (theta^(-2)) *(kappa- rho*theta*u*1i -d).*(1-exp(-d*T))./(1-g.*exp(-d*T)));
end
