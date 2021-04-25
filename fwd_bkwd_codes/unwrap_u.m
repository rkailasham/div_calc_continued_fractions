function [u_numeric] = unwrap_u(u_hat,transp_flag,L,varphi,ndim,Q,normQ,N)

%Unwraps the matrix "u" and returns its numeric value
%for a given chain configuration. The matrix "v" is obtained
%as an output of either "u_hat_jl_calc" or "transpose_u_hat_lj"
%"transp_flag" indicates if the input tensor is "u_hat" or
%its transpose.
%If this flag is 0, it means "u_hat" has been input.
%If this flag is 1, it means "u_hat_transpose" has been input.

%%% Here is the structure of "u" which has 10 columns
%%% Column 1  : ret_val from defn of U, rouse_prefac
%%% Column 2  : ret_val from legality of \bm{Q}_l
%%% Column 3  : iv_exp, exponent to iv prefactor
%%% Column 4  : seq_type (1) for "I", (2) for "D"
%%% Column 5  : Numerator subscript
%%% Column 6  : Denominator subscript
%%% Column 7  : lsq_start, starting index for continued product of L_i's
%%% Column 8  : lsq_fin, finishing index for continued product of L_i's
%%% Column 9  : first index, "f", and 
%%% Column 10 : second index,"s"
%%%             of the dyad (\bm{Q}_f\bm{Q}_s)/(Q_fQ_s)


p=(varphi/((2*varphi)+1))^2;
iv_pf=(varphi/((2*varphi)+1));
u_numeric=zeros(ndim,ndim);

[nrow,ncol]=size(u_hat);

for r=1:nrow
    u_scalar_numeric=unwrap_scalar_u_row(u_hat(r,:),...
        transp_flag,L,p,iv_pf,N);
    f=u_hat(r,(ncol-1));
    s=u_hat(r,(ncol));
    dyad=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));
    term=u_scalar_numeric*dyad;
    u_numeric=u_numeric+term;
end


end

