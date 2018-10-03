function h = Payoff_discounted(r,gamma,T,S,K,M,N)
    dt=T/N;
    h=zeros(N,M); 
    h(N,:)=exp(-r*T)*max(S(N+1,:)-K*exp(gamma*T),0);
    for i=1:N
        n=i+1;
        time=(n-1)*dt; 
        h(i,:)=exp(-r*time)*max(S(n,:)-K*exp(gamma*time),0);
    end
end