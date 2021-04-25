function [ret_val] = ret_rouse_val(i,j)
%Compares "i" and "j" and returns 
%--> "2" if i==j
%--> "-1" if |i-j|==1
%--> "0" for all other "i" and "j"

if(i==j)
    ret_val=2;
elseif(abs(i-j)==1)
    ret_val=-1;
else
    ret_val=0;
end

end

