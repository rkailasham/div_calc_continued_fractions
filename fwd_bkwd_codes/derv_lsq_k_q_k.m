function [partial] = derv_lsq_k_q_k(L,Q,normQ,k,N)
%Calculates the partial derivative of L^2_{k} w.r.t \bm{Q}_k
if(k==N)
    partial=0;
else
    pf=(2.*L(k))/(normQ(k+1)*normQ(k)*normQ(k));
    t1=normQ(k)*Q(k+1,:);
    t2=L(k)*normQ(k+1)*Q(k,:);
    partial=pf*(t1-t2);
end

end

