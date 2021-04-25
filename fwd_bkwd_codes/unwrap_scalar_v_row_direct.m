function [v_scalar_numeric] = unwrap_scalar_v_row_direct(v_hat_rownum,transp_flag,...
    L,M,P,iv_pf)

%Unwraps the scalar part of a given row (rownum) of matrix "v" and 
%returns its numeric value
%for a given chain configuration. 
%The suffix "direct" indicates that all quantities are evaluated
%in terms of "M" and "P", without recourse to the polynomial 
%expressions.
%The matrix "v" is obtained
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



[nrow,ncol]=size(v_hat_rownum);

if(~transp_flag)
    kcol=ncol-1;
else
    kcol=ncol;
end

multi_pf=v_hat_rownum(1)*v_hat_rownum(2)*v_hat_rownum(3);
ivfac=(iv_pf)^(v_hat_rownum(4));

if(v_hat_rownum(5)==1)
    ratio_seq=cont_prod_inv_series((v_hat_rownum(6)+1),(v_hat_rownum(7)),M);
elseif(v_hat_rownum(5)==2)
    ratio_seq=cont_prod_inv_series((v_hat_rownum(7)),(v_hat_rownum(6)-1),P);
end
lseq=cont_prod(v_hat_rownum(8),v_hat_rownum(9),L);
k=v_hat_rownum(kcol);
pf=1./(1.-M(k)-P(k));

v_scalar_numeric=multi_pf*pf*ivfac*ratio_seq*lseq;

end

