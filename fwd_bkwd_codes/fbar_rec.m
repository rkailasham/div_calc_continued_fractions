function [val] = fbar_rec(m,l,j,ell,compare_array,L)

if (m==1)
    val=1;
else
    val=0;
    for s=(j+2):l
        if(~ismember(s,compare_array))
            %ismember checks if "s" is present in
            %"compare_array" or not.
            val=val+(L(s)*L(s)*fbar_rec(m-1,l+2,s,ell,compare_array,L));
        end
    end
end

end

