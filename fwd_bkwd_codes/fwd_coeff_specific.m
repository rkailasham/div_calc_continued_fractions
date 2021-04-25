function [ret_val] = fwd_coeff_specific(p,L,n)
%returns the forward substitution coefficients
M=zeros(n,1);
if (n==1)
    M(n)=0;
else
    for i=2:n
        M(i)=(p*L(i-1)*L(i-1))/(1-M(i-1));
    end
end
ret_val=M(n);
end

