clc
clear
close all

% Defining Fuel as C_m H_n
m = 11;
n = 22;

% Calculating fuel mass flows

m_dot_air = [14.3 29.5 37.8];       % kg/s
P3 = 1e3 * [950 1950 3400];         % Pa
T3 = [620 765 850];                 % K
T4 = [1350 1580 1820];              % K

c_p3 = [1059 1094 1113];            % J/kg/K
c_p4 = [1191 1221 1244];            % J/kg/K

c_pavg = (c_p3 + c_p4) / 2;         % J/kg/K

h = 43e6;                           % J/kg
comb_efficiency = 1;                % -

m_dot_fuel = c_pavg .* m_dot_air .* (T4 - T3) ./ ( h * comb_efficiency - c_pavg .* T4 );


% Calculating equivalence ratio
% m_dot_air = 0.8 * m_dot_air;

m_dot_oxygen = 0.21 * m_dot_air;
fuel_to_oxg =  m_dot_fuel ./ m_dot_oxygen;
fuel_to_oxg_st = 2*(m * 12 + n * 1)/(33 * 2 * 16);

equivalence_ratio = fuel_to_oxg ./ fuel_to_oxg_st;
AFR = m_dot_air ./ m_dot_fuel;

fprintf('------------------------------\n')
for i = 1:3
    fprintf('Point %.0f: \n', i)
    fprintf('------------------------------\n')
    fprintf('Fuel Mass Flow: %f \n', m_dot_fuel(i))
    fprintf('Equivalence Ratio: %f \n', equivalence_ratio(i))
    fprintf('AFR: %f \n', AFR(i))
    fprintf('------------------------------\n')
end