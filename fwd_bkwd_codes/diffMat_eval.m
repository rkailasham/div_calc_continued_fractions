function [diffMat] = diffMat_eval(varphi,L,Q,normQ,N,ndim)
%Constructing the diffusion block matrix, \bm{\mathcal{D}}
ndiff=ndim*N;
diffMat=zeros(ndiff,ndiff);
fac=(varphi/((2*varphi)+1));

for j=1:N
    for k=1:N
        [v_hat] = v_hat_jl_calc(j,k,N);
        [v_numeric] = unwrap_v(v_hat,0,L,varphi,ndim,Q,normQ,N);
        idmat=ret_rouse_val(j,k)*eye(ndim); %identity matrix
        a=idmat-(fac*v_numeric);
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        diffMat(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=a;
    end
end

end

