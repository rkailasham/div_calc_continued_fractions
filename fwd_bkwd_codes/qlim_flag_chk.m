function [ret_val] = qlim_flag_chk(n,N)

%Returns "1" if the input index "n" is a part
%of the chain. Returns "0" otherwise"

if(n>=1&&n<=N)
    ret_val=1;
else
    ret_val=0;
end

end

