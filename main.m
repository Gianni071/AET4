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
X3 = 0.15;

rich.m_dot_fuel = m_dot_fuel;                                   % All incoming fuel
rich.m_dot_air = (0.2 + 2*X1) * m_dot_air;                      % Amount of air coming in
rich.AFR = rich.m_dot_air ./ rich.m_dot_fuel;
rich.phi = AFR_st ./ rich.AFR;

quench.m_dot_fuel = m_dot_fuel - rich.m_dot_air/AFR_st;         % All remaining fuel after the rich burn
quench.m_dot_air = (2*X2 + 0.2) * m_dot_air;                            % All oxygen has been used in rich so only 2*X2 air is coming in
quench.AFR = quench.m_dot_air ./ quench.m_dot_fuel;
quench.phi = AFR_st ./ quench.AFR;

lean.m_dot_fuel = quench.m_dot_fuel;                            % Same as quench since no fuel was used there
lean.m_dot_air = quench.m_dot_air + 2*X3*m_dot_air;   % All remaining air is burned now
lean.AFR = lean.m_dot_air ./ lean.m_dot_fuel;
lean.phi = AFR_st ./ lean.AFR;

%% Calculate Tad for all points

% Point 1
for i = drange(1:1:3)
    rich.T_ad(i) = FlameTemp(T3(i),rich.phi(i));
    quench.T_ad(i) = FlameTemp(T3(i),quench.phi(i));
    lean.T_ad(i) = FlameTemp(T3(1),lean.phi(i));
end

%% Residence time

R = 287;

rich.m_dot = rich.m_dot_air + rich.m_dot_fuel;
quench.m_dot = quench.m_dot_air + quench.m_dot_fuel;
lean.m_dot = lean.m_dot_air + lean.m_dot_fuel;

t_res = residence_time(rich.m_dot, rich.T_ad, P3) + residence_time(quench.m_dot, quench.T_ad, P3) + residence_time(lean.m_dot, lean.T_ad, P3);

heat_of_reaction = 6715.75;                                                 % kJ/mole fuel
M = 0.15429;                                                                % kg/mole
can_volume = 0.012;                                                         % m^3
heat_density_abs = m_dot_fuel ./ M .* heat_of_reaction ./ can_volume;       % kW/m^3
heat_density_norm = heat_density_abs ./ P3 .* 1e5;                          % kW/m^3/bar






fprintf('------------------------------\n')
for i = 1:3
    fprintf('Point %.0f: \n', i)
    fprintf('------------------------------\n')
    fprintf('Fuel Mass Flow: %f kg/s \n', m_dot_fuel(i))
    fprintf('Overall equivalence Ratio: %f \n', phi_overall(i))
    fprintf('AFR: %f \n', AFR(i))
    fprintf('------------------------------\n')
    fprintf('Equivalence ratio rich zone: %f \n', rich.phi(i))
    fprintf('Equivalence ratio quench zone: %f \n', quench.phi(i))
    fprintf('Equivalence ratio lean zone: %f \n', lean.phi(i))
    fprintf('------------------------------\n')
    fprintf('Rich adiabatic flame temp: %f K \n', rich.T_ad(i))
    fprintf('Quench adiabatic flame temp: %f K \n', quench.T_ad(i))
    fprintf('Lean adiabatic flame temp: %f K \n', lean.T_ad(i))
    fprintf('------------------------------\n')
    fprintf('Residence time: %f ms \n', 1000 .* t_res(i))
    fprintf('Absolute Heat Density: %f kW/m^3 \n', heat_density_abs(i))
    fprintf('Normalized Heat Density: %f kW/m^3/bar \n', heat_density_norm(i))
    fprintf('\n')
end

%% Plot Flame Temp as a function of equivalence ratio (T3 = 298)

Tarray = [];
phiarray = [];

for j = drange(0:0.01:2)
    Tad = FlameTemp(765,j);
    Tarray = [Tarray Tad];
    phiarray = [phiarray j];
end


plot(phiarray,Tarray)

