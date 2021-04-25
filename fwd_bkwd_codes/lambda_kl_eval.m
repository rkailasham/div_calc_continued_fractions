function [lambda_kl] = lambda_kl_eval(k,l,L,varphi,ndim,Q,normQ,N)
%Returns Lambda^{(k)}_l

lambda_kl=zeros(1,ndim);

p=(varphi/((2*varphi)+1))^2;
iv_pf=(varphi/((2*varphi)+1));
M_k = fwd_coeff_specific(p,L,k);
P_k = bkwd_coeff_specific(p,L,k,N);
pf=(1./(1.-M_k-P_k));

%By construction k,l go from 1 to N
%I have not implemented safety checks on 
%the values of k and l, but the function that 
%calls lambda_kl ought to be prudent.

if(l<k)
    [ret_val,iv_exp,seq_type,num_I,denom_I,...
        lsq_start,lsq_fin,vec_index] = gamma_kl(k,l,N);
    
    if(ret_val==1)
        ivfac=(iv_pf)^(iv_exp);        
        if(seq_type==1)
            %This check might feel redundant here,
            %but will be useful later on, while unwrapping
            ratio_i=fwd_coeff_poly(num_I,L,p,N)/fwd_coeff_poly(denom_I,L,p,N);
        else
            ratio_i=0;
        end
        lseq=cont_prod(lsq_start,lsq_fin,L);
        unit_vec=Q(vec_index,:)./normQ(vec_index);
        lambda_kl=ivfac*ratio_i*lseq*unit_vec;
    else
        lambda_kl(1,:)=0.; 
    end
    
elseif(l==k)
    unit_vec=Q(l,:)./normQ(l);
    lambda_kl=unit_vec;

elseif(l>k)    
    [ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = rho_kplus1_l((k+1),l,N);
    
    if(ret_val==1)
        ivfac=(iv_pf)^(iv_exp);
        if(seq_type==2)
            %This check might feel redundant here,
            %but will be useful later on, while unwrapping
            ratio_d=bkwd_coeff_poly(num_D,L,p,N)/bkwd_coeff_poly(denom_D,L,p,N);
        else
            ratio_d=0;
        end
        lseq=cont_prod(lsq_start,lsq_fin,L);
        unit_vec=Q(vec_index,:)./normQ(vec_index);
        lambda_kl=ivfac*ratio_d*lseq*unit_vec;
    else
        lambda_kl(1,:)=0.; 
    end
    
end
 
lambda_kl=pf*lambda_kl;

end

