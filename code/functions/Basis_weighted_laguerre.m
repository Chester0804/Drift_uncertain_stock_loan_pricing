% J: number of basis functions
function Y=Basis_weighted_laguerre(X,J)
    Y=zeros(J,length(X));
    for j=1:J 
        Y(j,:)=exp(-X./2).*Laguerre(X,j-1);
    end
end
               

                                              