function [ret_val,iv_exp,seq_type,num_I,...
    denom_I,lsq_start,lsq_fin,vec_index] = gamma_kl(k,l,N)
%Returns the parameters needed to construct Gamma^{(k)}_l
iv_exp=(k-l);
%seq_type=1 indicates that the sequence is 
%constructed using "I"
seq_type=1;
num_I=(l-1);
denom_I=(k-1);
lsq_start=l;
lsq_fin=(k-1);
vec_index=l;
flag_i_lim=i_n_compare(num_I,N);
flag_q_lim=qlim_flag_chk(vec_index,N);

if(flag_i_lim&&flag_q_lim&&(l<=k)&&(l>=1))
    ret_val=1;
else
    ret_val=0;
end

end

