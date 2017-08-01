clc;
% Stratification Monte Carlo Simulation
disp('Stratified Sampling')
n = 1000; n_st = 10; p = 1/n_st; m = 0; sq = 0;
for i = 0:n_st-1
    X = (rand(n,1)+i)*p;
    F = exp(X.^2);
    m = m + p*mean(F);
    sq = sq + std(F)^2*p^2/n;
end
disp('Est_Strat = '), disp(m)
disp('Var_Strat = '), disp(sq)

% Pure Monte Carlo Simulation
disp('Naive Monte-Carlo')
X=rand(n*n_st,1);
F=exp(X.^2);
m = mean(F);
sq = std(F)^2/(n*n_st);
disp('Est_Naive = '), disp(m)
disp('Var_Naive = '), disp(sq)