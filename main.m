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

AFR_st = 33 * 2 * 16 / ( 2*(m * 12 + n * 1) ) / 0.23;

AFR = m_dot_air ./ m_dot_fuel;
phi_overall = AFR_st ./ AFR;




% Calculating Rich, Quench and Lean zones equivalence ratios
X1 = 0.03;
X2 = 0.12;

rich.m_dot_fuel = m_dot_fuel;                                   % All incoming fuel
rich.m_dot_air = (0.2 + 2*X1) * m_dot_air;                      % Amount of air coming in
rich.AFR = rich.m_dot_air ./ rich.m_dot_fuel;
rich.phi = AFR_st ./ rich.AFR;

quench.m_dot_fuel = m_dot_fuel - rich.m_dot_air/AFR_st;         % All remaining fuel from the rich burn
quench.m_dot_air = 2*X2 * m_dot_air;                            % All oxygen has been used in rich so only 2*X2 air is coming in
quench.AFR = quench.m_dot_air ./ quench.m_dot_fuel;
quench.phi = AFR_st ./ quench.AFR;

lean.m_dot_fuel = quench.m_dot_fuel;                            % Same as quench since no fuel was used there
lean.m_dot_air = (1-0.2)*m_dot_air - rich.m_dot_air - quench.m_dot_air;   % All remaining air is burned now
lean.AFR = lean.m_dot_air ./ lean.m_dot_fuel;
lean.phi = AFR_st ./ lean.AFR;



fprintf('------------------------------\n')
for i = 1:3
    fprintf('Point %.0f: \n', i)
    fprintf('------------------------------\n')
    fprintf('Fuel Mass Flow: %f \n', m_dot_fuel(i))
    fprintf('Overall equivalence Ratio: %f \n', phi_overall(i))
    fprintf('AFR: %f \n', AFR(i))
    fprintf('------------------------------\n')
    fprintf('Equivalence ratio rich zone: %f \n', rich.phi(i)')
    fprintf('Equivalence ratio quench zone: %f \n', quench.phi(i)')
    fprintf('Equivalence ratio lean zone: %f \n', lean.phi(i)')
    fprintf('------------------------------\n')
    fprintf('\n')
end