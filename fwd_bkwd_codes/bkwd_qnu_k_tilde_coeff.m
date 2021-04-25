function [bkwd_qnu_k_tilde] = bkwd_qnu_k_tilde_coeff(nu,compare_array,k,L,p,minval,N)

if(nu<k||nu>=N)
    bkwd_qnu_k_tilde=0;
else
    bkwd_qnu_k_tilde=0;
    n=ceil((N-k)/2);
    fin=N+1;
    for mu=2:n
        init=k+(2*(mu-3));
        sum_term=gtilde_rec(mu-1,init,fin,nu,compare_array,L,minval);
        sgn=((-1)^mu);
        bkwd_qnu_k_tilde=bkwd_qnu_k_tilde+(sgn*((p^mu)*sum_term));
    end
end

end

