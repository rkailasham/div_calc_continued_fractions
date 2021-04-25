function [v_hat_jl_t] = transpose_v_hat_jl(v_hat_jl)

%Re-arranges the information stored in "v_hat_jl", so that its
%transpose may be assembled.
%This would involve exchanging the 10th and 11th columns of
% "v_hat_jl", so that the tensor (\bm{Q}_f\bm{Q}_s)/(Q_fQ_s) becomes
%(\bm{Q}_s\bm{Q}_f)/(Q_sQ_f). 

%Once the exchange of the two columns
%has been completed, remember that the prefactor,
%(1/(1-M(k)-P(k)) would depend now on the contents of the 11th column.

%"v_hat_jl" has 11 columns

v_hat_jl_t=v_hat_jl;
[nrow,ncol]=size(v_hat_jl);
temp1=v_hat_jl(:,ncol-1);
temp2=v_hat_jl(:,ncol);

v_hat_jl_t(:,ncol-1)=temp2;
v_hat_jl_t(:,ncol)=temp1;


end

