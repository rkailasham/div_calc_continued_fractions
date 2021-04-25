function [compare_array] = gen_list(varargin)
%Takes in an arbitrary number of arguments.
%For every arguement "ell" in the list
%returns an array containing [(ell-1) (ell) (ell+1)].
%Removes duplicate entries for simplicity

compare_array=[];
ell_list=cell2mat(varargin);

for i=1:nargin
    compare_array=[compare_array,(ell_list(i)-1),(ell_list(i)),(ell_list(i)+1)];
end

compare_array=unique(compare_array);
    
end

