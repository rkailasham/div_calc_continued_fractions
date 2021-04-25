function [ret_val] = bkwd_coeff_specific(p,L,n,N)
%returns the nth backward substitution coefficient
P = bkwd_coeff_all(p,L,N);
ret_val=P(n);
end

