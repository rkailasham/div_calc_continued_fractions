function [qnu_k_tilde] = bkwd_coeff_lsquared_nu_composite(nu,compare_array,k,L,p,minval,N)
%This returns the coefficient of L_{nu}^2 in D_k.
%These coefficients follow the composite rule.
%One can build up D_k by multiplying these
%coefficients with qnu_k_tilde and summing up.

if(nu<k||nu>=N)
    qnu_k_tilde=0;
else
    qnu_k_tilde=0;
    n=ceil((N-k)/2);
    fin=N+1;
    for mu=1:n
        init=k+(2*(mu-2));
        sum_term=gtilde_rec(mu,init,fin,nu,compare_array,L,minval);
        sgn=((-1)^mu);
        qnu_k_tilde=qnu_k_tilde+(sgn*((p^mu)*sum_term));
    end
end

end

