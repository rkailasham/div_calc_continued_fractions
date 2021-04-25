function [lambda_hat_kl] = lambda_hat_kl_calc(k,l,N)

%Returns the information necessary for constructing 
%the vector \hat{\Lambda}^{(k)}_l. There are 3 branches, depending
%on if l<k,l==k or l>k. The vector to be returned will 
%have 9 columns ---> 8 columns from the calculation of appropriate
%branch + 1 dummy column which will indicate the subscript of the
%vector which will premultiply \hat{\Lambda}^{(k)}_l.

lambda_hat_kl=zeros(1,9);

if(l<k)
    [ret_val,iv_exp,seq_type,num_I,denom_I,...
        lsq_start,lsq_fin,vec_index] = gamma_kl(k,l,N);
    lambda_hat_kl(1:7)=[ret_val,iv_exp,seq_type,num_I,denom_I,...
        lsq_start,lsq_fin];
    lambda_hat_kl(9)=vec_index;
    
elseif(l==k)
    lambda_hat_kl(1)=1;     %ret_val
    lambda_hat_kl(2)=0;     %iv_exp
    lambda_hat_kl(3)=1;     %seq_type (doesn't matter because ratio is 1)
    lambda_hat_kl(4)=1;     %num_I
    lambda_hat_kl(5)=1;     %denom_I
    lambda_hat_kl(6)=1;     %lsq_start
    lambda_hat_kl(7)=0;     %lsq_fin
    lambda_hat_kl(9)=l;     %vec_index

elseif(l>k)    
    [ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = rho_kplus1_l((k+1),l,N);
    lambda_hat_kl(1:7)=[ret_val,iv_exp,seq_type,num_D,denom_D,...
        lsq_start,lsq_fin];
    lambda_hat_kl(9)=vec_index;
    
end

lambda_hat_kl(8)=0; %This column will get filled with the index of the premultiplying 
                    %\bm{Q}_k. Use will become apparent when constructing \alpha.

end

