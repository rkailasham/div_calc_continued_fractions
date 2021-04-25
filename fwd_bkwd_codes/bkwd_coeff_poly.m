function [D_k] = bkwd_coeff_poly(k,L,p,N)
if(k>=1&&k<=(N+1))
    D_k=1;
else
    D_k=0;
    return;
end
n=ceil((N-k)/2);
fin=N-1;
for i=1:n
    init=k+(2*(i-1));
    sum_term=sum_call_g(i,init,fin,L);
    sgn=((-1)^i);
    D_k=D_k+(sgn*((p^i)*sum_term));
end

end

