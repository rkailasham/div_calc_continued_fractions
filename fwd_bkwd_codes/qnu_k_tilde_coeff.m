function [qnu_k_tilde] = qnu_k_tilde_coeff(nu,compare_array,k,L,p,minval)

if(nu<1||nu>=k)
    qnu_k_tilde=0;
else
    qnu_k_tilde=0;
    n=floor(k/2);
    init=-1;
    for mu=2:n
        fin=k-((2*mu)-5);
        sum_term=ftilde_rec(mu-1,fin,init,nu,compare_array,L,minval);
        sgn=((-1)^(mu));
        qnu_k_tilde=qnu_k_tilde+(sgn*((p^(mu))*sum_term));
    end    
end

end

