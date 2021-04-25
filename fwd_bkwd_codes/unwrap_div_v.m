function [div_v_numeric] = unwrap_div_v(v_hat,j,transp_flag,L,varphi,ndim,Q,normQ,N)

%Takes "v_hat" or its transpose as input, unwraps it,
%and returns its divergence w.r.t \bm{Q}_j as the numeric output for 
%for a given chain configuration. 
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
div_v_numeric=zeros(1,ndim);

[nrow,ncol]=size(v_hat);

for r=1:nrow
    [div_v_row_numeric] = unwrap_div_v_rowwise(v_hat(r,:),transp_flag,...
    j,L,Q,normQ,p,iv_pf,N,ndim);
    div_v_numeric=div_v_numeric+div_v_row_numeric;
end


end

