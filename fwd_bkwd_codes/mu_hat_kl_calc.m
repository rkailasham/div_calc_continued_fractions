function [mu_hat_kl] = mu_hat_kl_calc(k,l,N)
%Assembles the information necessary for
%constructing the tensor, mu_hat_kl.
%Multiplying this tensor by (1/1-M(k)-P(k)),
%a scalar, results in mu_kl.

vec_pre_mult=zeros(4,1);
vec_pre_mult(:,:)=k;
%This is the subscript for the vector \bm{Q}_k
%which pre-multiplies j_hat_kl.
%This information is stored as the second-to-last column
%in mu_hat_kl;
j_hat_kl = j_hat_kl_calc(k,l,N);
mu_hat_kl=j_hat_kl;
mu_hat_kl(:,9)=vec_pre_mult;

%remember that the size of mu_hat_kl is (4,10)
end

