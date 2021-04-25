function [ret_val] = i_n_compare(n,N)
%sets the return value based on 
%the subscript "n" of I_{n}.

if(n>=0&&n<=N)
    ret_val=1;
else
    ret_val=0;
    return;
end

end

