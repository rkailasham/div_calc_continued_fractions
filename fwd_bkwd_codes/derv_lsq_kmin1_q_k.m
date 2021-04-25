function [partial] = derv_lsq_kmin1_q_k(L,Q,normQ,k,N)
%Calculates the partial derivative of L^2_{k-1} w.r.t \bm{Q}_k
if(k==1)
    partial=0;
else
    pf=(2.*L(k-1))/(normQ(k-1)*normQ(k)*normQ(k));
    t1=normQ(k)*Q(k-1,:);
    t2=L(k-1)*normQ(k-1)*Q(k,:);
    partial=pf*(t1-t2);
end

end


