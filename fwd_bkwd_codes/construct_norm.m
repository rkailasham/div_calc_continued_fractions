function [normQ] = construct_norm(Q,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
normQ=zeros(N,1);
for i=1:N
    normQ(i)=norm(Q(i,:));
end

end

