function r = residence_time(m_dot_air, T_adiabatic_flame, P)
    R = 287;
    V = 0.004;
    r = P .* V ./(m_dot_air .* R .* T_adiabatic_flame);
end