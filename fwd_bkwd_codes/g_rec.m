function [val] = g_rec(m,l,j,L)

if (m==1)
    val=1;
else
    val=0;
    for s=l:(j-2)
        val=val+(L(s)*L(s)*g_rec(m-1,l-2,s,L));
    end
end

end

