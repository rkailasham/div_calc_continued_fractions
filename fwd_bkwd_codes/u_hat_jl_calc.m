function [u_hat_jl] = u_hat_jl_calc(j,l,N)

%Returns the information needed for assembling the tensor 
%U_hat_{jl}=2\alpha_hat_{jl}-\alpha_hat_{j-1,l}-\alpha_hat_{j+1,l}
%I persist with including the phrase "hat" here, because
%the premultiplication with (1/(1-M(k)-P(k)) will only be 
%performed during the unwrapping stage.

pf1=2*ones(1,1);
pf2=-1*ones(1,1);
pf3=pf2;

term1 = [pf1,alpha_hat_kl_calc(j,l,N)];
term2 = [pf2,alpha_hat_kl_calc((j-1),l,N)];
term3 = [pf3,alpha_hat_kl_calc((j+1),l,N)];

u_hat_jl=[term1;term2;term3];
[nrow,ncol]=size(u_hat_jl);

zero_ret_val_count=0;
for row=1:nrow
    chk=qlim_flag_chk(u_hat_jl(row,(ncol-1)),N);
    if(~chk||(u_hat_jl(row,2)==0))
        u_hat_jl(row,2)=0; %set ret_val to 0 if pre-multiplying
                           %vector is not legal.
        zero_ret_val_count=zero_ret_val_count+1;                   
    end
end     

%Lastly, it makes sense to only store columns with a non-zero ret_value

tempmat=u_hat_jl;
u_hat_jl=zeros((nrow-zero_ret_val_count),ncol);
build_count=0;

for row=1:nrow
    if(tempmat(row,2))
        build_count=build_count+1;
        u_hat_jl(build_count,:)=tempmat(row,:);
    end
end

%Typical size of "u_hat_jl" is 3x10, i.e; 3 rows and 10 columns

end

