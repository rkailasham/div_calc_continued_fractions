function [inv_pi_mk] = calc_inv_pi_mk(i,k,L,p,N)
%Returns the value of [1-M(i)-P(k)]^{-1}
t1=fwd_coeff_poly(i,L,p,N)/fwd_coeff_poly(i-1,L,p,N);
t2=bkwd_coeff_poly(k,L,p,N)/bkwd_coeff_poly(k+1,L,p,N);
sum=t1+t2-1.;
inv_pi_mk=1./sum;
end

