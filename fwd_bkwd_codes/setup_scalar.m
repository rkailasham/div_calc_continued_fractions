function [term] = setup_scalar(i,k,j,m,n,N,L,p)
M_i = fwd_coeff_specific(p,L,i);
P_k = bkwd_coeff_specific(p,L,k,N);
term=(1./(1.-M_i-P_k))*L(j)*L(j-1)*...
    (bkwd_coeff_poly(m,L,p,N)/bkwd_coeff_poly(n,L,p,N));    
end

