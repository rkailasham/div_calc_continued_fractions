function [L] = constructL(Q,normQ,N)
%Construct "L" which is the cosine of the angle
%between adjacent connector vectors
L=zeros(N-1,1);

for i=1:(N-1)
    cur=Q(i,:);
    nxt=Q(i+1,:);
    L(i)=(dot(cur,nxt))/(normQ(i)*normQ(i+1));
end

end

