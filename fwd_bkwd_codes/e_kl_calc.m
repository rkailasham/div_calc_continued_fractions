function [e_kl] = e_kl_calc(k,l,N)

%Returns the information necessary for constructing the vector,
%E^{(k)}_l=2\Gamma^{(k)}_{l}-\Gamma^{(k)}_{l-1}-\Gamma^{(k)}_{l+1}
%3 rows to indicate that E^{(k)}_l could have at most 3 terms.
%10 columns --> accomodate 1 prefactor + 8 return values from
%              gamma_kl(k,l,N) + one dummy column which will
%              indicate the subscript of the pre-multiplying \bm{Q}

e_kl=zeros(3,10);
rouse_prefac=[2;-1;-1;];
e_kl(:,1)=rouse_prefac;

[ret_val,iv_exp,seq_exp,num_I,...
denom_I,lsq_start,lsq_fin,vec_index] = gamma_kl(k,l,N);
e_kl(1,2:8)=[ret_val,iv_exp,seq_exp,num_I,denom_I,lsq_start,lsq_fin];
e_kl(1,10)=vec_index;

[ret_val,iv_exp,seq_exp,num_I,...
denom_I,lsq_start,lsq_fin,vec_index] = gamma_kl(k,(l-1),N);
e_kl(2,2:8)=[ret_val,iv_exp,seq_exp,num_I,denom_I,lsq_start,lsq_fin];
e_kl(2,10)=vec_index;

[ret_val,iv_exp,seq_exp,num_I,...
denom_I,lsq_start,lsq_fin,vec_index] = gamma_kl(k,(l+1),N);
e_kl(3,2:8)=[ret_val,iv_exp,seq_exp,num_I,denom_I,lsq_start,lsq_fin];
e_kl(3,10)=vec_index;

e_kl(:,9)=0; %This column will get filled with the index of the premultiplying 
             %\bm{Q}_k. Use will become apparent when constructing \mu.

end

