function [u_jl] = u_jl_eval(j,l,L,varphi,ndim,Q,normQ,N)

%Returns the tensor U_{jl}=2\alpha_{jl}-\alpha_{j-1,l}-\alpha_{j+1,l}
%Important to check that every element on the RHS exists.
%e.g U_{j-1,l} has no meaning when j=1;
%and U_{j+1,l} has no meaning when j=N;
%This check is implemented in the call to "alpha_kl_calc"

u_jl=2*alpha_kl_eval(j,l,L,varphi,ndim,Q,normQ,N)...
     -alpha_kl_eval((j-1),l,L,varphi,ndim,Q,normQ,N)...
     -alpha_kl_eval((j+1),l,L,varphi,ndim,Q,normQ,N);
    
end

