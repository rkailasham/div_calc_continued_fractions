function [T] = outerprod(v1,v2,n)
%returns the outer or dyadic product of two vectors "v1" and "v2",
%each of size "n".

T=zeros(n,n);

for i=1:n
    for j=1:n
        T(i,j)=v1(i)*v2(j);
    end
end

end

