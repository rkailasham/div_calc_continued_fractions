function [I_k] = fwd_coeff_poly(k,L,p,N)
if(k>=0&&k<=N)
    I_k=1;
else
    I_k=0;
    return;
end
n=floor(k/2);
init=1;
for i=1:n
    fin=k-((2*i)-1);
    sum_term=sum_call_f(i,init,fin,L);
    sgn=((-1)^i);
    I_k=I_k+(sgn*((p^i)*sum_term));
end

end

