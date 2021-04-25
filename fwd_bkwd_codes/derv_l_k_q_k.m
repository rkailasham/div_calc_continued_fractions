function [partial] = derv_l_k_q_k(L,Q,normQ,k,N)
%Calculates the partial derivative of L_{k} w.r.t \bm{Q}_k
if(k==N)
    partial=0;
else
    pf=(1./(normQ(k)));
    t1=(Q(k+1,:)./normQ(k+1));
    t2=L(k)*(Q(k,:)./normQ(k));
    partial=pf*(t1-t2);
end

end

