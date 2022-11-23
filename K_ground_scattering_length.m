%% Load Constants
cm = 1e-2;                  % centimeter [m]
cL = 299792458;             % speed of light [m/s]
c = 299792458;              % speed of light [m/s]
A0 = 1e-10;                 % Angstrom [m]

me = 9.1093837015e-31;      % electron mass [kg]
hbar = 1.054571817e-34;     % reduced Planck's constant [J s]
h = hbar*(2*pi);            % Planck's constant [J s]
a0 = 5.29177210903e-11;     % Borh radius in [m]
Eh = hbar^2/(me*a0^2);      % Hartree energy [J]
e0 = 4.3597447222071e-18;   % Hartree (energy) in J

amu = 1.66054e-27;          % atomic mass units [k]


m39 = 38.963706487*amu;     % 39K mass [kg]
m40 = 39.9637*amu;          % 40K mass [kg]

% Choose the mass
% m = m39;    
m = m40;    

% Reduced mass
mu = m*m/(2*m);

if m==m39
    mass_str = '$39\mathrm{K}$';
else
    mass_str = '$40\mathrm{K}$';    
end
    
%% Load Ground State Potentials Parameters
% Potentials are specified with position of angstroms.

% Spatial vector
RR = linspace(1,10000,2e5)';

% Potential for X^1 sigma^+_g (ground state) in m^-1
K_X_1_sigma;
funcs;
U_X_1_sigma = U(RR*A0);

% Potential for a^3 sigma^+_u (first excited) in m^-1
K_a_3_sigma;
funcs;
U_a_3_sigma = U(RR*A0);

% Normal energy with respect to the angstrom
U_xsigma_normalized = (U_X_1_sigma*cL*h)/(hbar^2/(2*mu*A0^2));
U_asigma_normalized = (U_a_3_sigma*cL*h)/(hbar^2/(2*mu*A0^2));


%% 
l=0;

Escat = 1e-3;

kVec = linspace(1e-3,.015,100)';

dVec = zeros(length(eVec),1);
for kk=1:length(eVec)
    disp(kk);
    E = kVec(kk)^2;
%     delta = solveScatteringPhase(RR,U_xsigma_normalized,E,l);
        delta = solveScatteringPhase(RR,U_asigma_normalized,E,l);

    dVec(kk) = delta;
end

%%


figure(4)
clf
a=54;
% a=70;
% p = linspace(

kt=linspace(.0009,0.015,100);
plot(kt,cotd,'r-','linewidth',2);
hold on
plot(kVec,cot(dVec),'ko')
hold on

re=40;
% re=35;
cotd=-1./(kt*a) + 0.5*re*kt;


xlabel('momentum');
ylabel('cot(delta)');