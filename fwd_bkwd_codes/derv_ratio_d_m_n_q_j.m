function [partial] = derv_ratio_d_m_n_q_j(m,n,j,N,L,p,Q,normQ)
%returns the derivative of (D_m/D_n) w.r.t \bm{Q}_j

D_m=bkwd_coeff_poly(m,L,p,N);
D_n=bkwd_coeff_poly(n,L,p,N);

derv_D_m=derv_d_k_q_j(m,j,L,p,Q,normQ,N);
derv_D_n=derv_d_k_q_j(n,j,L,p,Q,normQ,N);

t1=D_n*derv_D_m;
t2=D_m*derv_D_n;
denom=D_n*D_n;

partial=(t1-t2)/denom;

end

