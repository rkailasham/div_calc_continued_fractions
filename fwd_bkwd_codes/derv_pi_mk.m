function [partial,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N)
%Calculates the derivative of [1/(1-M_i-P_k)] w.r.t \bm{Q}_j

t1=fwd_coeff_poly(i,L,p,N)/fwd_coeff_poly(i-1,L,p,N);
derv_t1=derv_ratio_i_m_n_q_j(i,i-1,j,N,L,p,Q,normQ);

t2=bkwd_coeff_poly(k,L,p,N)/bkwd_coeff_poly(k+1,L,p,N);
derv_t2=derv_ratio_d_m_n_q_j(k,k+1,j,N,L,p,Q,normQ);

denom=(t1+t2-1);
num=-(derv_t1+derv_t2);

partial=num/(denom*denom);
end

