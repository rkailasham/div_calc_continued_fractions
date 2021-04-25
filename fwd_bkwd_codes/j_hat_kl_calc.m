function [j_hat_kl] = j_hat_kl_calc(k,l,N)
%Constructs \hat{J}^{(k})_l. This has four branches

term1=zeros(3,10);
term2=zeros(1,10);
term_prefac=(1-kron_delta(k,N));
kplus1=k+1;
%10 columns --> accomodate 1 prefactor + 9 return values so 
%              that the overall structure is similar to e_kl  or g_kl.

if((l<k)&&(l>=1))
    term1=e_kl_calc(k,l,N);
    term2(:,:)=0;
elseif(l==k)
    term1=e_kl_calc(k,l,N);
    term2(1,1)=-term_prefac;
    term2(1,2)=term_prefac;  %ret_val
    term2(1,3)=1;            %iv_exp
    term2(1,4)=2;            %seq_type
    term2(1,5)=(k+2);        %num_D
    term2(1,6)=(k+1);        %denom_D
    term2(1,7)=k;            %lsq_start     
    term2(1,8)=k;            %lsq_fin
    term2(1,9)=0;            %dummy; useful when building \mu
    term2(1,10)=(k+1);       %vec_index
elseif(l==(k+1))
    term1=g_kplus1l_calc(kplus1,l,N);
    term1(:,1)=term_prefac*term1(:,1);
    term2(1,1)=-term_prefac;
    term2(1,2)=term_prefac;  %ret_val
    term2(1,3)=0;            %iv_exp
    term2(1,4)=2;            %seq_type
    term2(1,5)=(N+1);        %num_D
    term2(1,6)=(N);          %denom_D
    term2(1,7)=1;            %lsq_start     
    term2(1,8)=0;            %lsq_fin
    term2(1,9)=0;            %dummy; useful when building \mu
    term2(1,10)=(k);         %vec_index  
    %Remember: there are no "L" terms in term2
    %which is why lsq_fin<lsq_start by construction.
elseif((l>(k+1))&&(l<=N))
    term1=g_kplus1l_calc(kplus1,l,N);
    term2(:,:)=0;    
end

j_hat_kl=[term1;term2];
end
