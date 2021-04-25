function [ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = y_kplus1_l(kplus1,m,N)
%Returns the parameters needed to construct Y^{(k+1)}_m
k=kplus1-1;
iv_exp=(m);
%seq_type=2 indicates that the sequence is 
%constructed using "D"
seq_type=2;
num_D=(m+k+1);
denom_D=(k+1);

%The starting and end-points of the continued
%product in "L" have been adjusted, by subtracting 1
%from both these limits. This means that "cont_prod"
%can be invoked by simply providing these modifed limits and
%"L" as input.

lsq_start=k;
lsq_fin=(k+m-1);

vec_index=(m+k);
flag_d_lim=d_n_compare(num_D,N);
flag_q_lim=qlim_flag_chk(vec_index,N);

if(flag_d_lim&&flag_q_lim&&(m>=1))
    ret_val=1;
else
    ret_val=0;
end

end

