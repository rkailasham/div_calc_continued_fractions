function [partial] = derv_l_kmin1_l_k_q_k(L,Q,normQ,k,N)
%Calculates the partial derivative of [L_{k-1}L_{k}] w.r.t \bm{Q}_k
if(k==1||k==N)
    partial=0;
else
    t1=normQ(k-1)*normQ(k)*L(k-1)*Q(k+1,:);
    t2=2*normQ(k-1)*normQ(k+1)*L(k-1)*L(k)*Q(k,:);
    t3=normQ(k)*normQ(k+1)*L(k)*Q(k-1,:);
    denom=normQ(k-1)*normQ(k)*normQ(k)*normQ(k+1);
    partial=(t1-t2+t3)/denom;
end

end


