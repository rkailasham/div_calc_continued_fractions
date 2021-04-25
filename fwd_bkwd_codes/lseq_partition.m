function [pf2,term3,derv3] = lseq_partition(j,lsq_start,lsq_fin,L,Q,normQ,N)
%Function to determine which term(s) in the continued product 
%are to be taken out of the derivative, and which term(s)
%are to be kept inside the bracket. "j" is the subscript in 
%\bm{Q}_{j} with respect to which the derivative is taken.

%This function has been built to cover the most general case.
%In actual use, it is possible that several of these branches 
%may not be triggered.

nterms=(lsq_fin-lsq_start)+1;
lseq=cont_prod(lsq_start,lsq_fin,L);

%%% In every branch, note that the following conservation rule
%%% applies : pf2*term3=lseq

if(nterms<=0)
    term3=lseq;
    derv3=0;    
    
elseif(nterms==1)
    if(lsq_start==(j-1))
        term3=lseq; 
        derv3=derv_l_kmin1_q_k(L,Q,normQ,j,N);
    elseif(lsq_start==j)
        term3=lseq;
        derv3=derv_l_k_q_k(L,Q,normQ,j,N);
    else
        term3=1;
        derv3=0;
    end
    
elseif(nterms>1)
    
    if((j>lsq_start)&&(j<=(lsq_fin)))%product of 2 terms in "term3"
        term3=L(j-1)*L(j);
        derv3=derv_l_kmin1_l_k_q_k(L,Q,normQ,j,N);
    elseif(j==lsq_start)%only one term in "term3", effectively
        term3=L(lsq_start);
        derv3=derv_l_k_q_k(L,Q,normQ,j,N);
    elseif(j==(lsq_fin+1))%only one term in "term3", effectively      
        term3=L(lsq_fin);
        derv3=derv_l_kmin1_q_k(L,Q,normQ,j,N);
    else
        term3=1;
        derv3=0;
    end

end

pf2=lseq/term3;

end

