function [] = Vega_Asian_Pathwise_v_LikelihoodRatio()
% This code implements the pathwise and likelihood ratio methods for estimating the vega of an Asian option
% in the GBM Black-Scholes framework. We plot the variance of the 2 estimators as a function of the number
% of weekly observations used to compute the average. See Figure 7.3 in Glasserman's Monte-Carlo book.
% E.G 7.3.3 on pg 405 contains LR vega score and e.g. 7.2.3 on pg 390 contains the pathwise estimator.

% Input Parameters:  X0 = initial stock price, T = maturity, K = strike, sig = vol, r = riskfree rate = drift
S0 = 100;       K = 100;      r = .05;    sig = .3;
%randn('state',sum(100*clock));
%rng('default') ; % Sets random number generator seed to default value

Maturities = 1:1:256;  %2.^[1:8];    % Maturities of Asian options
NumAsians = length(Maturities);  % No. of Asian options
n = max(Maturities);     % Longest maturity (in weeks)
dt = 1/52;  % Corresponding to weekly averaging
T = Maturities*dt;   % Vector of maturities in years. (Assuming each option expires in last averaging week)
NumSamples = 500000;  % Number of sample paths. Should take around 1 million
Con1 = (r - sig^2/2)*dt;   Con2 = sig*sqrt(dt);
S = zeros(NumSamples,n);     dS_dSig = zeros(NumSamples,n);   ScoreTerms = zeros(NumSamples,n);
SumZ = zeros(NumSamples,1);

for i=1:n
    Z = randn(NumSamples,1);
    if (i==1)
        S(:,i) = S0 * exp(Con1 + Con2*Z);
    else
        S(:,i) = S(:,i-1) .*exp(Con1 + Con2*Z);
    end
    % Next two terms are requried for the pathwise vega calculation
    SumZ = SumZ + sqrt(dt)*Z;
    dS_dSig(:,i) = S(:,i).* (-sig*Maturities(i)*dt + SumZ);
    
    % Next term is requried for the LR vega calculation
    ScoreTerms(:,i) = (Z.^2 - 1)/sig - sqrt(dt)*Z;
end

% Now compute the pathwise vega's mean and variance
CumSum_dS_dSig = cumsum(dS_dSig, 2);
AsianAverage = cumsum(S,2) ./ repmat(Maturities,NumSamples,1);
PathwiseVegaSamples = CumSum_dS_dSig .* [AsianAverage > K] .* repmat(exp(-r*T) ./ Maturities, NumSamples, 1);
PathwiseVega_Mean = mean(PathwiseVegaSamples);
PathwiseVega_Var = var(PathwiseVegaSamples);   

% Now compute the LR vega's mean and variance
AsianPayoffDiscounted = max(0, AsianAverage - K).* repmat(exp(-r*T), NumSamples, 1);
Score = cumsum(ScoreTerms,2);
LRVegaSamples = AsianPayoffDiscounted .* Score;
LRVega_Mean = mean(LRVegaSamples);
LRVega_Var = var(LRVegaSamples);

% Now plot the estimated vega for these Asian options for each method as a function of weeks
plot(Maturities,[PathwiseVega_Mean; LRVega_Mean],'linewidth',1);
hold on
grid on
set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Weeks','fontsize',15)
ylabel('Estimated Vega','rot',90,'fontsize',15)
h_legend = legend('Pathwise','Likelihood Ratio');
set(h_legend,'FontSize',13);

% Now plot the variance (per replication) of the vega for each method as a function of weeks
figure
semilogy(Maturities,[PathwiseVega_Var; LRVega_Var],'linewidth',1);
hold on
grid on
set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Weeks','fontsize',15)
ylabel('Variance','rot',90,'fontsize',15)
h_legend = legend('Pathwise','Likelihood Ratio');
set(h_legend,'FontSize',13);

%%
% Now plot mean +- 1.96 \hat{sigma} / sqrt(NumSamples) for each method
CI_LR = [LRVega_Mean - 1.96*sqrt(LRVega_Var)/sqrt(NumSamples);  LRVega_Mean; LRVega_Mean + 1.96*sqrt(LRVega_Var)/sqrt(NumSamples)];  
figure 
plot(Maturities,CI_LR,'linewidth',1);
hold on
grid on
set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Weeks','fontsize',15)
ylabel('Estimated Vega','rot',90,'fontsize',15)
title('Estimated Vega with 95% confidence bounds via the LR Method')

CI_Pathwise = [PathwiseVega_Mean - 1.96*sqrt(PathwiseVega_Var)/sqrt(NumSamples);  PathwiseVega_Mean; PathwiseVega_Mean + 1.96*sqrt(PathwiseVega_Var)/sqrt(NumSamples)];  
figure 
plot(Maturities,CI_Pathwise,'linewidth',1);
hold on
grid on
set(gca, 'linewidth',1,'fontsize',10); 
xlabel('Number of Weeks','fontsize',15)
ylabel('Estimated Vega','rot',90,'fontsize',15)
title('Estimated Vega with 95% confidence bounds via the Pathwise Method')
end




