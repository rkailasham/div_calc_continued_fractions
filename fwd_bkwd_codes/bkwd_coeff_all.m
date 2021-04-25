function [P] = bkwd_coeff_all(p,L,N)
%returns the backward substitution coefficients
P=zeros(N,1);
P(N)=0.;
for i=(N-1):-1:1
    P(i)=(p*L(i)*L(i))/(1.-P(i+1));
end

end

