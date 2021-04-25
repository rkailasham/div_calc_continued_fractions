function [qnu_k_tilde] = fwd_coeff_lsquared_nu_composite(nu,compare_array,k,L,p,minval)
%This returns the coefficient of L_{nu}^2 in I_k.
%These coefficients follow the composite rule.
%One can build up I_k by multiplying these
%coefficients with qnu_k_tilde and summing up.

if(nu<1||nu>=k)
    qnu_k_tilde=0;
else
    qnu_k_tilde=0;
    n=floor(k/2);
    init=-1;
    for mu=1:n
        fin=k-((2*mu)-3);
        sum_term=ftilde_rec(mu,fin,init,nu,compare_array,L,minval);
        sgn=((-1)^mu);
        qnu_k_tilde=qnu_k_tilde+(sgn*((p^mu)*sum_term));
    end
end

end

