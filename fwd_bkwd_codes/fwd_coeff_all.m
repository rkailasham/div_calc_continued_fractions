function [M] = fwd_coeff_all(p,L,N)
%returns the forward substitution coefficients
M=zeros(N,1);
M(1)=0.;
for i=2:N
    M(i)=(p*L(i-1)*L(i-1))/(1.-M(i-1));
end

end

