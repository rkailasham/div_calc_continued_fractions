function [cp] = cont_prod_inv_series(start,fin,A)
%returns the continued product of (1-A_i)^{-1} multiplied
%from i=start to i=fin
cp=1;

for i=start:fin
    val=1./(1.-A(i));
    cp=cp*val;
end

end

