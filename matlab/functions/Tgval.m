function Tg = Tgval(T1g, Lu)
% TGVAL - Tg = Tgval(T1g, Lu)
%   evaluates point values of T[g_m](nu) at N equispaced points where
%   T[g_m] is defined by
%   
%     T^{(1)}[g](nu) + T^{(1)}_0[g]*Lu(nu)
%     
    n = length(T1g);
    L0 = -1/n*sum( Lu );
    T10g = 1/(n*L0)*sum( T1g );
    Tg = T1g + T10g*Lu;
end
