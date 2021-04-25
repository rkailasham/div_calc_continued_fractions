function [term] = setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p)
dyad=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));
M_i = fwd_coeff_specific(p,L,i);
P_k = bkwd_coeff_specific(p,L,k,N);
term=(1./(1.-M_i-P_k))*L(j)*L(j-1)*...
    (bkwd_coeff_poly(m,L,p,N)/bkwd_coeff_poly(n,L,p,N))*dyad;    
end

