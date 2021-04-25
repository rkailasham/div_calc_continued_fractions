function [term] = sum_call_g(mu,init,fin,L)
term=0;

for i=init:fin
    term=term+(L(i)*L(i)*g_rec(mu,(init-2),i,L));
end

end

