function[Tad] = FlameTemp(T3,phi)

%% Reactants enthalpies from T3 (Interpolation)
kcaltokj = 4.184;
dHf1 = 25.19*kcaltokj; %kJ/mol
dHf2 = 40.55*kcaltokj; %kJ/mol
dHf3 = 50.39*kcaltokj; %kJ/mol
dOx1 = 9.892; %kJ/mol
dOx2 = 14.67; %kJ/mol
dOx3 = 17.54; %kJ/mol
Hffuel298 = -65.114*kcaltokj;
R = 8.3144598;

%% Oxygen polynomial coefficients
Aoxarray = [31.32234	30.03235	20.91111];
Boxarray = [-20.23531	8.772972	10.72071];
Coxarray = [57.86644	-3.988133	-2.020498];
Doxarray = [-36.50624	0.788313	0.146449];
Eoxarray = [-0.007374	-0.741599	9.245722];
Foxarray = [-8.903471	-11.32468	5.337651];
Goxarray = [246.7945	236.1663	237.6185];

%% CO2 polynomial coefficients
Acoarray = [24.99735	58.16639];
Bcoarray = [55.18696	2.720074];
Ccoarray = [-33.69137	-0.492289];
Dcoarray = [7.948387	0.038844];
Ecoarray = [-0.136638	-6.447293];
Fcoarray = [-403.6075	-425.9186];
Gcoarray = [228.2431	263.6125];
Hcoarray = [-393.5224	-393.5224];

%% H2O polynomial coefficients
Awatarray = [30.09200	41.96426];
Bwatarray =	[6.832514	8.622053];
Cwatarray =	[6.793435	-1.499780];
Dwatarray =	[-2.534480	0.098119];
Ewatarray =	[0.082139	-11.15764];
Fwatarray =	[-250.8810	-272.1797];
Gwatarray =	[223.3967	219.7809];
Hwatarray =	[-241.8264	-241.8264];

%% N2 polynomial coefficients
An2array =[ 28.98641	19.50583	35.51872];
Bn2array =	[1.853978	19.88705	1.128728];
Cn2array =	[-9.647459	-8.598535	-0.196103];
Dn2array =	[16.63537	1.369784	0.014662];
En2array =	[0.000117	0.527601	-4.553760];
Fn2array =	[-8.671914	-4.935202	-18.97091];
Gn2array =	[226.4168	212.3900	224.9810];

%% Fuel polynomial coefficients
Afarray = [0.47049127e+01 0.25897423e+02];
Bfarray = [0.58911327e-01  0.55462092e-01];
Cfarray = [0.10500014e-03 -0.17337738e-04];
Dfarray = [-0.18501088e-06 0.17582452e-08];
Efarray = [0.82238431e-10 0.63971899e-13];
Ffarray = [-0.37391863e+05 -0.46337805e+05];
Gfarray = [0.13439449e+02 -0.11004780e+03];

%% Calculations
heat_of_reaction = 6715.75; %kJ/mol fuel

if T3 == 620
    dHf = dHf1;
    dOx = dOx1;
end

if T3 == 765
    dHf = dHf2;
    dOx = dOx2;
end

if T3 == 850
    dHf = dHf3;
    dOx = dOx3;
end

Tarr = [];
diffarr = [];
for T = drange(298:1:3000)
    % N2 pols
    if T < 2000
        N2pols = [An2array(2) Bn2array(2) Cn2array(2) Dn2array(2) En2array(2) Fn2array(2) Gn2array(2)];
    end
    if T > 2000
        N2pols = [An2array(3) Bn2array(3) Cn2array(3) Dn2array(3) En2array(3) Fn2array(3) Gn2array(3)];
    end
    
    % H2O pols
    if T < 1700
        Watpols = [Awatarray(1) Bwatarray(1) Cwatarray(1) Dwatarray(1) Ewatarray(1) Fwatarray(1) Gwatarray(1) Hwatarray(1)];
    end
    if T > 1700
        Watpols = [Awatarray(2) Bwatarray(2) Cwatarray(2) Dwatarray(2) Ewatarray(2) Fwatarray(2) Gwatarray(2) Hwatarray(2)];
    end
    
    % CO2 pols
    if T < 1200
        COpols = [Acoarray(1) Bcoarray(1) Ccoarray(1) Dcoarray(1) Ecoarray(1) Fcoarray(1) Gcoarray(1) Hcoarray(1)]; 
    end
    
    if T > 1200
        COpols = [Acoarray(2) Bcoarray(2) Ccoarray(2) Dcoarray(2) Ecoarray(2) Fcoarray(2) Gcoarray(2) Hcoarray(2)]; 
    end
    
    % O2 pols
    if T < 700
        Oxpols = [Aoxarray(1) Boxarray(1) Coxarray(1) Doxarray(1) Eoxarray(1) Foxarray(1) Goxarray(1)];
    end
    if T > 700 & T < 2000
        Oxpols = [Aoxarray(2) Boxarray(2) Coxarray(2) Doxarray(2) Eoxarray(2) Foxarray(2) Goxarray(2)];
    end
    if T > 2000
        Oxpols = [Aoxarray(3) Boxarray(3) Coxarray(3) Doxarray(3) Eoxarray(3) Foxarray(3) Goxarray(3)];
    end
    
    % Fuel pols
    if T < 1000
        Fuelpols = [Afarray(1) Bfarray(1) Cfarray(1) Dfarray(1) Efarray(1) Ffarray(1) Gfarray(1)];
    end
    
    if T > 1000
        Fuelpols = [Afarray(2) Bfarray(2) Cfarray(2) Dfarray(2) Efarray(2) Ffarray(2) Gfarray(2)];
    end
    % t for polynomials
    t = T/1000;
    
    % Number of moles
    nO2tot = 16.5/phi;
    nCO2 = 11/phi;
    nH2O = 11/phi;
    nN2 = 61.51/phi;
    nfuel = 1-(1/phi);
    
    % Heat needed to heat up products, leftover O2 and N2
    dHO2 = Oxpols(1)*t + Oxpols(2)*(t^2)/2 + Oxpols(3)*(t^3)/3 + Oxpols(4)*(t^4)/4 - Oxpols(5)/t + Oxpols(6);
    dHCO = COpols(1)*t + COpols(2)*(t^2)/2 + COpols(3)*(t^3)/3 + COpols(4)*(t^4)/4 - COpols(5)/t + COpols(6) - COpols(8);
    dHwat = Watpols(1)*t + Watpols(2)*(t^2)/2 + Watpols(3)*(t^3)/3 + Watpols(4)*(t^4)/4 - Watpols(5)/t + Watpols(6) - Watpols(8);
    dHN2 = N2pols(1)*t + N2pols(2)*(t^2)/2 + N2pols(3)*(t^3)/3 + N2pols(4)*(t^4)/4 - N2pols(5)/t + N2pols(6);
    dHfuel =  R*T*(Fuelpols(1) + Fuelpols(2)*T/2 + Fuelpols(3)*(T^2)/3 + Fuelpols(4)*(T^3)/4 + Fuelpols(5)*(T^4)/5 + Fuelpols(6)/T);
    dHfuel = dHfuel - (Hffuel298*1000);
    dHfuel = dHfuel/1000;
    
    dHN298 = 8.89;
    
    if phi > 1
        lefthand = heat_of_reaction/phi + dHf + nO2tot*dOx;
        righthand = nCO2*dHCO + nH2O*dHwat + nN2*dHN2 + nfuel*dHfuel;
    end
    
    if phi == 1
        lefthand = heat_of_reaction + dHf + nO2tot*dOx;
        righthand = nCO2*dHCO + nH2O*dHwat + nN2*dHN2;
    end
    
    if phi < 1
        nO2 = nO2tot - 16.5;
        lefthand = heat_of_reaction + dHf + nO2tot*dOx;
        righthand = nCO2*dHCO + nH2O*dHwat + nN2*dHN2 + nO2*dHO2;
    end
    
    diff = abs(righthand - lefthand);
    diffarr = [diffarr diff];
    Tarr = [Tarr T];
end
    
[discard,i] = min(diffarr);
Tad = Tarr(i);
