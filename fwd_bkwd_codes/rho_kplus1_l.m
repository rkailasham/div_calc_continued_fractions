function [ret_val,iv_exp,seq_type,num_D,...
    denom_D,lsq_start,lsq_fin,vec_index] = rho_kplus1_l(kplus1,l,N)
%Returns the parameters needed to construct rho^{(k+1)}_l
k=kplus1-1;
iv_exp=(l-k);
%seq_type=2 indicates that the sequence is 
%constructed using "D"
seq_type=2;
num_D=(l+1);
denom_D=(k+1);

%The starting and end-points of the continued
%product in "L" have been adjusted, by subtracting 1
%from both these limits. This means that "cont_prod"
%can be invoked by simply providing these modifed limits and
%"L" as input.

lsq_start=k;
lsq_fin=(l-1);

vec_index=l;
flag_d_lim=d_n_compare(num_D,N);
flag_q_lim=qlim_flag_chk(vec_index,N);

if(flag_d_lim&&flag_q_lim&&(l>=kplus1))
    ret_val=1;
else
    ret_val=0;
end

end

