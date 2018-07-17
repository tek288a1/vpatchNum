function [T1g, T10g] = T1gval(m, q0, omega0, dthdnu, A1g, c11g, Lu)
% T1GVAL - T1g = T1gval(m, omega0, dthdnu, A1g, c1g)
%   evaluates point values of T^{(1)}[g_m](nu) at N equispaced points
%   where T^{(1)}[g] is defined by
%   
%     q_0(nu)*( omega0(nu)*A_1[g](nu) - c_1^{(1)}[g]*dth/dnu(nu) )
%   
%   and the scaled average T^{(1)}_0[g_m] defined by
%   
%     1/(2*pi*L0) \int_{0}^{2*pi} T^{(1)}[g_m](nu) d nu
%
    n = length(omega0);
    T1g = q0.*( omega0.*A1g - c11g*dthdnu ); 
    L0 = -1/n*sum( Lu );
    T10g = 1/(n*L0)*sum( T1g );   
end

