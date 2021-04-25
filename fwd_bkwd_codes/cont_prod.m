function [cp] = cont_prod(start,fin,A)
%returns the continued product of A_i multiplied
%from i=start to i=fin
cp=1;

for i=start:fin
    val=A(i);
    cp=cp*val;
end

end

