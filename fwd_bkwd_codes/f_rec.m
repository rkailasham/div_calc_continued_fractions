function [val] = f_rec(m,l,j,L)

if (m==1)
    val=1;
else
    val=0;
    for s=(j+2):l
        val=val+(L(s)*L(s)*f_rec(m-1,l+2,s,L));
    end
end

end

