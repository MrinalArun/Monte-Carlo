clear;
clc;
lambda=5;
delta=1/100;
c=600;
A=10000;
N=100000;
alpha=0.95;

%solve for t*
t = delta - lambda/c;
%calculate M_t*
M_t=delta*lambda/((delta-t)*(lambda+t*c));
Z=zeros(N,1);
Count_step = floor(N/10);

p = (delta*lambda)/((lambda+delta*c)*M_t*(delta-t));
for i=1:N
    j=1;
    S=0;
    while S<A
        U=rand; %simulate X
        if U< p
            X=exprnd(1/(delta-t));
        else
            X=-exprnd(1/(lambda/c+t));
        end
        S=S+X;
        j=j+1;
    end
    Z(i)=exp(-t*S+j*log(M_t)); %estimator
    if mod(i,Count_step)==0 %track progress
        disp([num2str(i/N*100),' % completed.']);
    end
end

mean_Z=mean(Z);
disp('mean_z = '), disp(mean_Z)
std_Z=std(Z)/sqrt(N);
lb=mean_Z+norminv(1-alpha)*std_Z;
rb=mean_Z+norminv(alpha)*std_Z;
disp('CI = '), disp([lb, rb])