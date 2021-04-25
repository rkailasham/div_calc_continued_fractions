function [g_kplus1l] = g_kplus1l_calc(kplus1,l,N)

%Returns the information necessary for constructing the vector,
%G^{(k+1)}_l, where l=(k+1+m). The formula is
%G^{(k+1)}_l=2Y^{(k+1)}_{m+1}-Y^{(k+1)}_{m}-\Y^{(k+1)}_{m+2}
%3 rows to indicate that G^{(k+1)}_l could have at most 3 terms.
%10 columns --> accomodate 1 prefactor + 8 return values from
%              y_kplus1_l(kplus1,m,N) + one dummy column which will
%              indicate the subscript of the pre-multiplying \bm{Q}

k=kplus1-1;
m=l-k-1;

g_kplus1l=zeros(3,10);

rouse_prefac=[2;-1;-1;];
g_kplus1l(:,1)=rouse_prefac;

[ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = y_kplus1_l(kplus1,(m+1),N);
g_kplus1l(1,2:8)=[ret_val,iv_exp,seq_type,num_D,denom_D,lsq_start,lsq_fin];
g_kplus1l(1,10)=vec_index;

[ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = y_kplus1_l(kplus1,m,N);
g_kplus1l(2,2:8)=[ret_val,iv_exp,seq_type,num_D,denom_D,lsq_start,lsq_fin];
g_kplus1l(2,10)=vec_index;

[ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = y_kplus1_l(kplus1,(m+2),N);
g_kplus1l(3,2:8)=[ret_val,iv_exp,seq_type,num_D,denom_D,lsq_start,lsq_fin];
g_kplus1l(3,10)=vec_index;

g_kplus1l(:,9)=0; %This column will get filled with the index of the premultiplying 
                  %\bm{Q}_k. Use will become apparent when constructing \mu.

end

