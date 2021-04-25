function [v_numeric] = unwrap_v(v_hat,transp_flag,L,varphi,ndim,Q,normQ,N)

%Unwraps the matrix "v" and returns its numeric value
%for a given chain configuration. The matrix "v" is obtained
%as an output of either "v_hat_jl_calc" or "transpose_v_hat_lj"
%"transp_flag" indicates if the input tensor is "v_hat" or
%its transpose.
%If this flag is 0, it means "v_hat" has been input.
%If this flag is 1, it means "v_hat_transpose" has been input.

%%% Here is the structure of "v" which has 11 columns
%%% Column 1  : ret_val from defn of V, rouse_prefac
%%% Column 2  : ret_val from E or G, from j_hat_calc
%%% Column 3  : ret_val from legality of \bm{Q}_l
%%% Column 4  : iv_exp, exponent to iv prefactor
%%% Column 5  : seq_type (1) for "I", (2) for "D"
%%% Column 6  : Numerator subscript
%%% Column 7  : Denominator subscript
%%% Column 8  : lsq_start, starting index for continued product of L_i's
%%% Column 9  : lsq_fin, finishing index for continued product of L_i's
%%% Column 10 : first index, "f", and 
%%% Column 11 : second index,"s"
%%%             of the dyad (\bm{Q}_f\bm{Q}_s)/(Q_fQ_s)

p=(varphi/((2*varphi)+1))^2;
iv_pf=(varphi/((2*varphi)+1));
v_numeric=zeros(ndim,ndim);

[nrow,ncol]=size(v_hat);

for r=1:nrow
    v_scalar_numeric=unwrap_scalar_v_row(v_hat(r,:),...
        transp_flag,L,p,iv_pf,N);
    f=v_hat(r,(ncol-1));
    s=v_hat(r,(ncol));
    dyad=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));
    term=v_scalar_numeric*dyad;
    v_numeric=v_numeric+term;
end


end

