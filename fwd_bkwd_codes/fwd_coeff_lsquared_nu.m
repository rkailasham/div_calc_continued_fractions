function [qnu_k] = fwd_coeff_lsquared_nu(nu,compare_array,k,L,p)
%This returns the coefficient of L_{nu}^2 in I_k.
%However, these coefficients do not follow the composite rule.
%One cannot build up I_k by naively multiplying these
%coefficients with qnu_k and summing up: this would
%result in overcounting.

qnu_k=0;
n=floor(k/2);
init=-1;
for mu=1:n
    fin=k-((2*mu)-3);
    sum_term=fbar_rec(mu,fin,init,nu,compare_array,L);
    sgn=((-1)^mu);
    qnu_k=qnu_k+(sgn*((p^mu)*sum_term));
end

end

