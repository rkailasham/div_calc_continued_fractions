function [partial] = derv_ratio_i_m_n_q_j(m,n,j,N,L,p,Q,normQ)
%returns the derivative of (I_m/I_n) w.r.t \bm{Q}_j

I_m=fwd_coeff_poly(m,L,p,N);
I_n=fwd_coeff_poly(n,L,p,N);

derv_I_m=derv_i_k_q_j(m,j,L,p,Q,normQ,N);
derv_I_n=derv_i_k_q_j(n,j,L,p,Q,normQ,N);

t1=I_n*derv_I_m;
t2=I_m*derv_I_n;
denom=I_n*I_n;

partial=(t1-t2)/denom;

end

