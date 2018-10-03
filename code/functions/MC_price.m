% This is the function for pricing stock loan by Monte Carlo simulation
function SL_price=MC_price(r,gamma,a,b,Delta,S0,pi0,K,T)
N=2^6; %number of time step
M=1e4; %number of paths
J=5; %number of basis function
dt=T/N; %step size
t=(0:1:N)*dt;
pi=zeros(N+1,M);
S=zeros(N+1,M);

% Generate stock price paths
pi(1,:)=pi0 * ones(1,M);
for n=2:N+1
    dW=sqrt(dt)*randn(1,M);
    pi(n,:)=pi(n-1,:).*(1+(Delta*(1-pi(n-1,:))).*dW);
end
S(1,:)=S0 * ones(1,M);
for n=2:N+1
    dW=sqrt(dt)*randn(1,M);
    S(n,:)=S(n-1,:).*(1+dW+(Delta*pi(n-1,:)+b)*dt);
end

% Estimate beta
h=Payoff_discounted(r,gamma,T,S,K,M,N);
Beta = Regression_beta(S,h,N,M,J);

% Generate second set of stock price paths
pi_2(1,:)=pi0 * ones(1,M);
for n=2:N+1
    dW=sqrt(dt)*randn(1,M);
    pi_2(n,:)=pi_2(n-1,:).*(1+(Delta*(1-pi_2(n-1,:))).*dW);
end
S_2(1,:)=S0 * ones(1,M);
for n=2:N+1
    dW=sqrt(dt)*randn(1,M);
    S_2(n,:)=S_2(n-1,:).*(1+dW+(Delta*pi_2(n-1,:)+b)*dt);
end

% Calculate stock loan price
h_2=Payoff_discounted(r,gamma,T,S_2,K,M,N);
U=zeros(N,M);
c=zeros(N,M);
for i=1:N
    n=i+1;
    if (i==N)
        U(i,:)=h_2(i,:);
    else
        Psi=Basis_weighted_laguerre(S(n,:),J);
        c(i,:)=(Beta(:,i)')*Psi;
        U(i,:)=max(h(i,:),c(i,:));
    end
end
Vec_rule=zeros(1,M);
for m=1:M
    Vec_rule(m)=h(find(h(:,m)>=c(:,m),1),m);
end
SL_price=sum(Vec_rule)/M;
end