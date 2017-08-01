clc; clear;


A = [0.5 0 0; 0.3 0.6 0; 0.2 0.4 1];
B = [0.7 0.4 0.8; 0.3 0.6 0.2];
h0 = [0.9 0.1 0];
v = [1 2 1];

K = size(A,1);
T = size(v,2);

% alpha recursion

alpha = zeros(K,T);
for i=1:K
   alpha(i,1) = B(v(1),i)*h0(i); 
end

for j=2:T
   for i=1:K
      alpha(i,j) = B(v(j),i)*(A(i,:)*alpha(:,j-1)); 
   end
end

% beta recursion

beta = zeros(K,T); 
beta(:,T) = 1;

for j=T-1:-1:1
    for i=1:K
       beta(i,j) = B(v(j+1),:)*(A(:,i).*beta(:,j+1)); 
    end
end

display(['a) P(v(1:3)) = ',num2str(sum(alpha(:,T)))]);

display(['b) P(h(1)|v(1:3)) = [',num2str((alpha(:,1).*beta(:,1)/(alpha(:,1)'*beta(:,1)))'),']']);

% mu recursion

mu = zeros(K,T);
mu(:,T) = 1;

for j=T-1:-1:1
   for i=1:K
      mu(i,j) = max(B(v(j+1),:)'.*(A(:,i).*mu(:,j+1))); 
   end
end

hmax = zeros(1,T);
[~,hmax(1)] = max(B(v(1),:).*h0.*mu(:,1)');
for j=2:T
    [~,hmax(j)] = max(B(v(j),:).*A(:,hmax(j-1))'.*mu(:,1)');
end

display(['c) max. lik. path = [',num2str(hmax),']']);