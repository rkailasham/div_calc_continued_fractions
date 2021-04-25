function [alpha_kl] = alpha_kl_eval(k,l,L,varphi,ndim,Q,normQ,N)
%Returns the tensor \alpha_{kl}

flag_q_lim=qlim_flag_chk(k,N);

if(flag_q_lim==1)
    pf=Q(k,:)'./normQ(k);
    alpha_kl=pf*lambda_kl_eval(k,l,L,varphi,ndim,Q,normQ,N);
else
    alpha_kl=zeros(ndim,ndim);
end


end

