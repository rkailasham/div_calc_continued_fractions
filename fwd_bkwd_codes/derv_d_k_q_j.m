function [derv_d_k] = derv_d_k_q_j(k,j,L,p,Q,normQ,N)
%Calculates the derivative of D_k w.r.t the j-th
%connector vector
if(j<k)
    derv_d_k=0;
else
    nu=j-1;
    t2=bkwd_coeff_lsquared_nu_composite(nu,gen_list(nu),k,L,p,nu,N)*derv_lsq_kmin1_q_k(L,Q,normQ,j,N);
    nu=j;
    t3=bkwd_coeff_lsquared_nu_composite(nu,gen_list(nu),k,L,p,nu,N)*derv_lsq_k_q_k(L,Q,normQ,j,N);
    sum=0;
    for i=k:(j-2)
        term=L(i)*L(i)*derv_bkwd_qnu_k_tilde_q_j(i,j,k,L,p,Q,normQ,N);           
        sum=sum+term;
    end
    derv_d_k=sum+t2+t3;
end

end

