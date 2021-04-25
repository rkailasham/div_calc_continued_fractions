function [div_v_row_numeric] = unwrap_div_v_rowwise(v_hat_rownum,transp_flag,...
    j,L,Q,normQ,p,iv_pf,N,ndim)

%Unwraps the divergence of a given row (rownum) of matrix "v" and 
%returns its numeric value
%for a given chain configuration. 

%The formula is div(fT)=(grad(f)).T + f(div(T)) where "f"
%is scalar and "T" is a tensor

%This is evaluated by calculating the two terms on the RHS of
%the above equation separately, with "fac1=(grad(f)).T" and 
%"fac2=f(div(T))".


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

% % p=(varphi/((2*varphi)+1))^2;
% % iv_pf=(varphi/((2*varphi)+1));

div_v_row_numeric=zeros(1,ndim);

[nrow,ncol]=size(v_hat_rownum);

if(~transp_flag)
    kcol=ncol-1;
else
    kcol=ncol;
end

%Calculation of fac1=(grad(f)).T

multi_pf=v_hat_rownum(1)*v_hat_rownum(2)*v_hat_rownum(3);
ivfac=(iv_pf)^(v_hat_rownum(4));
pf1=multi_pf*ivfac;

k=v_hat_rownum(kcol);
[derv1,denom] = derv_pi_mk(k,k,j,p,L,Q,normQ,N);
piece1=(1./denom);

if(v_hat_rownum(5)==1)
    piece2=fwd_coeff_poly(v_hat_rownum(6),L,p,N)/fwd_coeff_poly(v_hat_rownum(7),L,p,N);
    derv2=derv_ratio_i_m_n_q_j(v_hat_rownum(6),v_hat_rownum(7),j,N,L,p,Q,normQ);
elseif(v_hat_rownum(5)==2)
    piece2=bkwd_coeff_poly(v_hat_rownum(6),L,p,N)/bkwd_coeff_poly(v_hat_rownum(7),L,p,N);
    derv2=derv_ratio_d_m_n_q_j(v_hat_rownum(6),v_hat_rownum(7),j,N,L,p,Q,normQ);
end

[pf2,piece3,derv3] = lseq_partition(j,v_hat_rownum(8),v_hat_rownum(9),L,Q,normQ,N);
prefac=pf1*pf2;
first=v_hat_rownum(10);
sec=v_hat_rownum(11);
dyad=(Q(first,:)'*Q(sec,:))./(normQ(first)*normQ(sec));

term1=piece1*piece2*derv3;
term2=piece1*piece3*derv2;
term3=piece2*piece3*derv1;

%(grad(f)).T

fac1=prefac*(term1+term2+term3)*dyad;

%Calculation of fac2=f(div(T))

f=unwrap_scalar_v_row(v_hat_rownum,transp_flag,L,p,iv_pf,N);
divT=derv_q_tens(first,sec,j,Q,normQ);

fac2=f*divT;

div_v_row_numeric=fac1+fac2;

end

