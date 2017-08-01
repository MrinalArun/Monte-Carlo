function [] = GibbsHierarchical()
% This function implements the posterior inference for the hierarchical model as outlined in section 11.6 of Gelman et al's 3rd
% edition of Bayesian Data Analysis
% See also my slides and the section on MCMC convergence diagnostics to see notation etc. We will use m/2 chains for a total of n0
% + 2n iterations. We will discard the first n0 samples from each chain and then split the remainder of each chain into chains so
% that after the burn-in we will have a total of 2m chains each of length n.

rng(20);
n0 = 100000;   n = 100000;
N = n0 + 2*n;
NumChains = 4;  % Number of chain that we run. Will split each into 2 chains after discarding burnin samples
m = 2*NumChains;

% Data
J = 4;  % # of groups
NumJ = [4,6,6,8]; 
nimax = max(NumJ);
Y = zeros(J,nimax);
Y(1,1:NumJ(1)) = [62, 60, 63, 59];
Y(2,1:NumJ(2)) = [63, 67, 71, 64, 65, 66];
Y(3,1:NumJ(3)) = [68, 66, 71, 67, 68, 68];
Y(4,1:NumJ(4)) = [56, 62, 60, 61, 63, 64, 63, 59];

% Initialize matrices that will hold samples from posterior
theta1 = zeros(n,m);    theta2 = zeros(n,m);    theta3 = zeros(n,m);    theta4 = zeros(n,m);
sigsq = zeros(n,m);     mu = zeros(n,m);    tausq = zeros(n,m);

% Now generate the samples, discard burnin samples and split remainder of each chain into 2
for i=1:NumChains  
    [theta_samps, mu_samps, sigsq_samps, tausq_samps] = GibbsSample(N,Y,J,NumJ);    
    
    theta1(:,2*i-1) = theta_samps(1,(n0+1):(n0+n))';    theta1(:,2*i) = theta_samps(1,(n0+n+1):N)';     
    theta2(:,2*i-1) = theta_samps(2,(n0+1):(n0+n))';    theta2(:,2*i) = theta_samps(2,(n0+n+1):N)'; 
    theta3(:,2*i-1) = theta_samps(3,(n0+1):(n0+n))';    theta3(:,2*i) = theta_samps(3,(n0+n+1):N)';
    theta4(:,2*i-1) = theta_samps(4,(n0+1):(n0+n))';    theta4(:,2*i) = theta_samps(4,(n0+n+1):N)';    
    sigsq(:,2*i-1) = sigsq_samps((n0+1):(n0+n))';       sigsq(:,2*i) = sigsq_samps((n0+n+1):N)';  
    mu(:,2*i-1) = mu_samps((n0+1):(n0+n))';             mu(:,2*i) = mu_samps((n0+n+1):N)';    
    tausq(:,2*i-1) = tausq_samps((n0+1):(n0+n))';       tausq(:,2*i) = tausq_samps((n0+n+1):N)';
end


%Now plot each posterior quantity
figure(1);      hold on;     subplot(2,2,1);    plot(theta1);       title('Theta1');
ylim([mean(theta1(:,1))-3*std(theta1(:,1)),mean(theta1(:,1))+3*std(theta1(:,1))]); 
subplot(2,2,2);     plot(theta2);       title('Theta2');
ylim([mean(theta2(:,1))-3*std(theta2(:,1)),mean(theta2(:,1))+3*std(theta2(:,1))]); 
subplot(2,2,3);     plot(theta3);       title('Theta3');
ylim([mean(theta3(:,1))-3*std(theta3(:,1)),mean(theta3(:,1))+3*std(theta3(:,1))]); 
subplot(2,2,4);     plot(theta4);       title('Theta4');
ylim([mean(theta4(:,1))-3*std(theta4(:,1)),mean(theta4(:,1))+3*std(theta4(:,1))]); 
hold off

figure(2);   hold on;    plot(mu);      title('mu');
ylim([mean(mu(:,1))-3*std(mu(:,1)),mean(mu(:,1))+3*std(mu(:,1))]); 
hold off

figure(3);  hold on;    plot(sigsq);    title('sigsq');
ylim([mean(sigsq(:,1))-3*std(sigsq(:,1)),mean(sigsq(:,1))+3*std(sigsq(:,1))]); 
hold off

figure(4);  hold on;    plot(tausq);    title('tausq');
ylim([mean(tausq(:,1))-3*std(tausq(:,1)),mean(tausq(:,1))+3*std(tausq(:,1))]); 
hold off

% Now perform Gelman-Rubin convergence diagnostics
R_theta1 = Gelman_Rubin_Diagnostics(theta1)
R_theta2 = Gelman_Rubin_Diagnostics(theta2)
R_theta3 = Gelman_Rubin_Diagnostics(theta3)
R_theta4 = Gelman_Rubin_Diagnostics(theta4)

R_sigsq = Gelman_Rubin_Diagnostics(sigsq)
R_mu = Gelman_Rubin_Diagnostics(mu)
R_tausq = Gelman_Rubin_Diagnostics(tausq)

% Now compute quantiles
p = [0.025 0.25 0.5 0.75 0.975];
theta1_q = quantile(theta1(:),p)
theta2_q = quantile(theta2(:),p)
theta3_q = quantile(theta3(:),p)
theta4_q = quantile(theta4(:),p)
mu_q = quantile(mu(:),p)
sig_q = quantile(sqrt(sigsq(:)),p)
tau_q = quantile(sqrt(tausq(:)),p)


end

%*************************************************************
function [theta_samps, mu_samps, sigsq_samps, tausq_samps] = GibbsSample(N,Y,J,NumJ)

% Initialize outputs
theta_samps = zeros(J,N);   mu_samps = zeros(1,N);      sigsq_samps = zeros(1,N);   tausq_samps = zeros(1,N); 
n = sum(NumJ);

% Generate widely dispersed estimates for the thetas and mu
for i=1:J
    RandInd = randi(NumJ(i),1); % Generates a random integer from 1:NumI(i)
    theta_samps(i,1) = Y(i,RandInd);  % This is indeed an overdispersed starting point for theta_i
end
mu_samps(1) = mean(theta_samps(:,1));


for k=2:N 
    % Generate tausq from Inv-Chi^2 (J-1, param_tau)
    param_tau = sum((theta_samps(:,k-1) - mu_samps(k-1)).^2) / (J-1);
    tausq_samps(k) = GenerateInvChi(J-1,param_tau);
    
    % Generate sigsq from Inv-Chi^2 (J-1n, param_sig)
    param_sig = 0;
    for j=1:J
       param_sig = param_sig + sum((Y(j,1:NumJ(j)) - theta_samps(j,k-1)).^2);
    end
    param_sig = param_sig/n;
    sigsq_samps(k) = GenerateInvChi(n,param_sig);
    
    % Generate mu
    mu_samps(k) = mean(theta_samps(:,k-1)) + sqrt(tausq_samps(k)/J) * randn;
    
    % Generate the theta_j's
    for j=1:J
        temp2 = 1/tausq_samps(k) + NumJ(j)/sigsq_samps(k);
        temp1 = mu_samps(k)/tausq_samps(k) + (NumJ(j)/sigsq_samps(k)) * mean(Y(j,1:NumJ(j)));
        param1 = temp1/temp2;
        param2 = 1/temp2;      
        theta_samps(j,k) = param1 + sqrt(param2)*randn;  
    end   
   
end

end

%*************************************************************
function [sample] = GenerateInvChi(p1,p2)

y = chi2rnd(p1);
sample = p1*p2/y;
end

%*************************************************************
function [R] = Gelman_Rubin_Diagnostics(PostSamples)
% PostSamples is an (n,m) matrix of posterior samples where m = # of chains and n = # samples in each chain

n = length(PostSamples(:,1));
avg = mean(PostSamples,1);
B = n*var(avg);
W = mean(var(PostSamples));
VarEst = (n-1)*W/n + B/n;
R = sqrt(VarEst/W);

end





    
    
