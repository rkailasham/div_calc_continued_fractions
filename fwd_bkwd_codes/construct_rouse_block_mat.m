function [rb_mat] = construct_rouse_block_mat(N,ndim)

ndsize=N*ndim;
rb_mat=zeros(ndsize,ndsize);

for j=1:N
    for k=1:N
        idmat=ret_rouse_val(j,k)*eye(ndim); %identity matrix
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        rb_mat(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=idmat;
    end
end

end

