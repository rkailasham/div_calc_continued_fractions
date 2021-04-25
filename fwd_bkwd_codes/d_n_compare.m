function [ret_val] = d_n_compare(n,N)
%sets the return value based on 
%the subscript "n" of D_n.

if(n>=1&&n<=(N+1))
    ret_val=1;
else
    ret_val=0;
    return;
end

end

