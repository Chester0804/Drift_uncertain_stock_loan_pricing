function Beta=Regression_beta(S,h,N,M,J)
    Beta=zeros(J,N);
    c=zeros(N,M);
    U=h(N,:);
    
    for i=(N-1):-1:1
        n=i+1;
        Psi=Basis_weighted_laguerre(S(n,:),J);
        
        sum1=0;
        sum2=0;
        
        for j=1:M
            sum1=sum1+Psi(:,j)*Psi(:,j)';
            sum2=sum2+Psi(:,j)*U(j);
        end
        
        Beta(:,i)=sum1\sum2;
        c(i,:)=(Beta(:,i)')*Psi;
        
        U=max(h(i,:),c(i,:));
        
    end
end    
        