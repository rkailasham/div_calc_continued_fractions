function [derv_qnu_k_tilde] = derv_bkwd_qnu_k_tilde_q_j(nu,j,k,L,p,Q,normQ,N)
%Calculates the derivative of the backward coefficient w.r.t the j-th
%connector vector
temp=j-nu;
if(j<(nu+2))
    derv_qnu_k_tilde=0;
else
    compare_array=gen_list(nu,j-1);
    minval=min(nu,j-1);
    if(temp>2)
        t1 = bkwd_qnu_k_tilde_coeff(nu,compare_array,k,L,p,minval,N)*derv_lsq_kmin1_q_k(L,Q,normQ,j,N);
    else
        t1=0;
    end
    compare_array=gen_list(nu,j);
    minval=min(nu,j);
    if(j==N)
        t2=0;
    else
        t2 = bkwd_qnu_k_tilde_coeff(nu,compare_array,k,L,p,minval,N)*derv_lsq_k_q_k(L,Q,normQ,j,N);
    end
    derv_qnu_k_tilde=t1+t2; 
end

end

