function dxdt = Zenteno_Model_2020_Full(t,x,p,int_time,exp_temp)

%Time varying temperature
T   = interp1(int_time,exp_temp,t)+273.15;     %Linear data interpolation

%State Variables

X    = x(1);           % Biomass (g/L)
N    = x(2);           % Nitrogen (g/L)
G    = x(3);           % Glucose (g/L)
F    = x(4);           % Fructose (g/L)
E    = x(5);           % Ethanol (g/L)
% CO2  = x(6);           % CO2 (g/L)
% Temp = x(7);           % must temeprature (K)

%Estimated Parameters
mu0    = p(1);        % 1/h (Boulton 1979)
betaG0 = p(2);        % kg E/kg bio/h (Salmon 2003)
betaF0 = p(3);        % kg E/kg bio/h (Salmon 2003)
Kn0    = p(4);        % kg N/m3 (Cramer 2002)
Kg0    = p(5);        % kg G/m3 (Cramer 2002)
Kf0    = p(6);        % kg F/m3 (Cramer 2002)
Kig0   = p(7);        % kg G/m3 (estimated salmon 2003, zenteno 2010)
Kie0   = p(8);        % kg E/m3 (Boulton 1979)
Kd0    = p(9);        % 1/h     (Zenteno 2010)

Yxn    = p(10);       % kg Bio/kg N (Zenteno 2010), biomass/nitrogen yield coeff
Yxg    = p(11);       % kg bio/kg G (Zenteno 2010), biomass/glucose yield coeff
Yxf    = p(12);       % kg bio/kg F (zenteno 2010), biomass/fructose yield coeff
Yeg    = p(13);       % kg E/ kg G (Zenteno 2010), ethanol/glucose yield coeff
Yef    = p(14);       % kg E/ kg F (Zenteno 2010), ethanol/fructose yield coeff

%Fixed Parameters
Cde     = 0.0415;     % m3/kg E (Salmon 2003)
Etd     = 130000.0;   % kj/kmol
V_must  = 0.003;      % m3
M_must  = V_must*1098;% kg
Fr      = 10;         % kg/h
Trin    = 15+273.15;
Tenv    = 30+273.15;  % environment temperature
%Afc     = 9.23;      % m2 heat transfer area of bioreactor/cooling system
Afc     = 0.094;      % m2 heat transfer area of bioreactor/cooling system
Ace     = 0.094;      % m2 heat transfer area of cooling system/environment
%Ufc     = 180;       % kJ/m2 h K, Global heat tranfer coeff bioreactor/cooling system
Ufc     = 16.72;      % kJ/m2 h K, Global heat tranfer coeff bioreactor/cooling system
Uce     = 16.72;      % kJ/m2 h K, Global heat tranfer coeff cooling system/environment
cpm     = 4.15;       % kJ/kg K (Boulton 1996) (average must and wine)
cpr     = 4.19;       % kJ/Kg K (Boulton et al. 1996)
R       = 8.314;      % kJ/kmol/K (Boulton 1979)
DHrx    = -651.28;    % kJ/kg    (Boulton 1979)
Eac     = 59453.0;    % kJ/kmol (Boulton 1979)
Eafe    = 11000.0;    % kJ/kmol (Zenteno 2010, trial and error)
EaKn    = 46055.0;    % kJ/kmol (Boulton 1979)
EaKg    = 46055.0;    % kJ/kmol (Boulton 1979)
EaKf    = 46055.0;    % kJ/kmol (Boulton 1979)
EaKig   = 46055.0;    % kJ/kmol (Boulton 1979)
EaKie   = 46055.0;    % kJ/kmol (Boulton 1979)
Eam     = 37681.0;    % kJ/kmol (Boulton 1979)
m0      = 0.01;       % kgS/kg bio/h
klatf   = 0.07;       % 1/h (Gee and ramirez 1994)
Csat    = 1.5;        % kg/m3 Saturation of dissolved CO2

%Constitutive Equations

% maximum specific growth rate, 1/h
mu_max    = mu0*exp(Eac*(T-300)/(300*R*T));
% maximum ethanol specific production rate from glucose, kg E/kg bio/h
betaG_max = betaG0*exp(Eafe*(T-296.15)/(296.15*R*T));
% maximum ethanol specific production rate from fructose, kg E/kg bio/h
betaF_max = betaF0*exp(Eafe*(T-296.15)/(296.15*R*T));
% Nitrogen compounds saturation constant, kg K/m3
Kn        = Kn0*exp(EaKn*(T-293.15)/(293.15*R*T));
% Glucose saturation constant, kg G/m3
Kg        = Kg0*exp(EaKg*(T-293.15)/(293.15*R*T));
% Fructose saturation constant, kg F/m3
Kf        = Kf0*exp(EaKf*(T-293.15)/(293.15*R*T));
% Glucose inhibition constant, kg G/m3
Kig       = Kig0*exp(EaKig*(T-293.15)/(293.15*R*T));
% Ethanol inhibition constant, kg E/m3
Kie       = Kie0*exp(EaKie*(T-293.15)/(293.15*R*T));
% maintenace coefficient, kg S/kg bio/h
m         = m0*exp(Eam*(T-293.3)/(293.3*R*T));
% growth rate, 1/h
mu        = mu_max*(N/(N+Kn));
% ethanol production rate from glucose
beta_G    = betaG_max*(G/(G+Kg))*(Kie/(E+Kie));
% ethanol production rate from fructose
beta_F    = betaF_max*(F/(F+Kf))*(Kig/(G+Kig))*(Kie/(E+Kie));
% threshold temperature of thermal death
Td        = -0.0001*E^3+0.0049*E^2-0.1279*E+315.89;

% specific death rate
if T>=Td
    Kd    = Kd0*exp((Cde*E)+(Etd*(T-305.65))/(305.65*R*T));
else
    Kd    = 0;
end

Ws        = (G/M_must)+(F/M_must);
We        = E/M_must;
fc        = -0.00075*T+1.2242;
rho_must  = (176.2*Ws^2-59.24*We^2 + 50.34*Ws*We+361.0*Ws-150.08*We + 998.69)*fc;
Tr        = (Fr*cpr*Trin+Ufc*Afc*T+Uce*Ace*Tenv)/((Fr*cpr)+(Ufc*Afc)+(Uce*Ace));


%Differential Equations
dXdt   = mu*X-Kd*X;                                 % Biomass
dNdt   = -mu*(X/Yxn);                               % Nitrogen
dGdt   = -((mu/Yxg)+(beta_G/Yeg)+m*(G/(G+F)))*X;    % Glucose
dFdt   = -((mu/Yxf)+(beta_F/Yef)+m*(F/(G+F)))*X;    % Fructose
dEdt   = (beta_G+beta_F)*X;                         % Ethanol

dxdt   = [dXdt;dNdt;dGdt;dFdt;dEdt];