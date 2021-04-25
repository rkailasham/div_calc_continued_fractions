function [val] = ftilde_rec(m,l,j,ell,compare_array,L,minval)

if (m==1)
    val=1;
else
    val=0;
    for s=(j+2):l
        if(~ismember(s,compare_array)&&s>minval)
            %ismember checks if "s" is present in
            %"compare_array" or not.
            val=val+(L(s)*L(s)*ftilde_rec(m-1,l+2,s,ell,compare_array,L,minval));
        end
    end
end

end

