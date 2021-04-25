function [qnu_k] = bkwd_coeff_lsquared_nu(nu,compare_array,k,L,p,N)
%This returns the coefficient of L_{nu}^2 in D_k.
%However, these coefficients do not follow the composite rule.
%One cannot build up D_k by naively multiplying these
%coefficients with qnu_k and summing up: this would
%result in overcounting.

qnu_k=0;
n=ceil((N-k)/2);
fin=N+1;
for mu=1:n
    init=k+(2*(mu-2));
    sum_term=gbar_rec(mu,init,fin,nu,compare_array,L);
    sgn=((-1)^mu);
    qnu_k=qnu_k+(sgn*((p^mu)*sum_term));
end

end

