function [derv_i_k] = derv_i_k_q_j(k,j,L,p,Q,normQ,N)
%Calculates the derivative of I_k w.r.t the j-th
%connector vector
if(j>k)
    derv_i_k=0;
else
    nu=j-1;
    t2=fwd_coeff_lsquared_nu_composite(nu,gen_list(nu),k,L,p,nu)*derv_lsq_kmin1_q_k(L,Q,normQ,j,N);
    nu=j;
    t3=fwd_coeff_lsquared_nu_composite(nu,gen_list(nu),k,L,p,nu)*derv_lsq_k_q_k(L,Q,normQ,j,N);
    sum=0;
    for i=1:(j-2)
        term=L(i)*L(i)*derv_qnu_k_tilde_q_j(i,j,k,L,p,Q,normQ,N);           
        sum=sum+term;
    end
    derv_i_k=sum+t2+t3;
end

end

