function [divQ] = derv_q_tens(i,k,j,Q,normQ)
%returns the divergence of [\bm{Q}_i\bm{Q}_k/(Q_iQ_k)]
%w.r.t \bm{Q}_j

if(j==i&&j~=k)
    divQ=(2./(normQ(k)*normQ(j)))*Q(k,:);
elseif(j==k&&j~=i)
    divQ=(2./(normQ(i)*normQ(j)))*Q(i,:);
elseif(j==i&&j==k)
    divQ=(2./(normQ(j)*normQ(j)))*Q(j,:);
else
    divQ=0;
end

end

