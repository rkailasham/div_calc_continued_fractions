function [term] = sum_call_f(mu,init,fin,L)
term=0;

for i=init:fin
    term=term+(L(i)*L(i)*f_rec(mu,(fin+2),i,L));
end

end

