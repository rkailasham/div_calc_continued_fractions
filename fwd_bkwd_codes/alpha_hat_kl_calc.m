function [alpha_hat_kl] = alpha_hat_kl_calc(k,l,N)

%Returns the information needed for constructing the 
%tensor \hat{\alpha}_{kl}. This has the same size as 
%lambda_hat_kl.

[lambda_hat_kl]=lambda_hat_kl_calc(k,l,N);
alpha_hat_kl=lambda_hat_kl;
alpha_hat_kl(8)=k;

end

