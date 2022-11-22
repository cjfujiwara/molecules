
clear all
%% Constants
a0=5.29177210903e-11; % bohr radius
A0 = 10^-10; % angstrom
cm=0.01;

%% Parameters

% Important separations
R0 = 2.55*A0;
R1 = 3.025*A0;
R2 = 18.5*A0;
R3 = 100*A0;

% Repulsive branch
A1 = 6.15621691e3*cm^-1;
C = 1.4954055;      % original value
C = 1.4944-0.00015945; % CF Values by tweaking short range

A2 = 5.7847022e4*cm^-1;

% Potential well
bVec=[1.110701931e4*cm^-1;
    -4.36925369044300638e1*cm^-1;
    4.76297891975524326e4*cm^-1;
    -2.09601115045173974e4*cm^-1;
    -7.01984407530179014e4*cm^-1;
    -3.01001775238950504e4*cm^-1;
    1.55968559608675539e6*cm^-1;
    -5.00914579944727710e5*cm^-1;
    -5.73511453432456404e7*cm^-1;
    1.15249130719885722e8*cm^-1;
    1.04092244644117713e9*cm^-1;
    -3.806592600279434456e9*cm^-1;
    -7.09216903641955090e9*cm^-1;
    5.45659790420072861e10*cm^-1;
    -5.92912726136713181e10*cm^-1;
    -2.36187587414028015e11*cm^-1;
    1.17116364854221021e12*cm^-1;
    -2.50404040001589111e12*cm^-1;
    -1.89576460672672681e12*cm^-1;
    2.86804447527311953e13*cm^-1;
    -6.17994630589855469e13*cm^-1;
    -3.03985720429939258e13*cm^-1;
    3.59974601606658937e14*cm^-1;
    -5.75771007718742750e14*cm^-1;
    -5.16723472499104531e13*cm^-1;
    1.47635288904728050e15*cm^-1;
    -2.40901143313356450e15*cm^-1;
    1.95020021872318475e15*cm^-1;
    -8.34049836318447375e14*cm^-1;
    1.50984509388737000e14*cm^-1];
Rm = 4.551*A0;
a=0.26;

% Long Range
C3s = 5.483104e5*cm^-1;
C6s = 5.122770e7*cm^-1;
C8s = 2.6654e9*cm^-1;

C3p = C3s/2;
C6p = 3.0319e7*cm^-1;
C8p = 1.02877e9*cm^-1;


D = 17474.5848*cm^-1;
nu = 13023.6587*cm^-1; % wavevector of 4s-4p
lambda_bar = 1/nu; % lambda bar is inverse wavevector?
Delta = 57.7103*cm^-1;

B2 = 6.0213760140;
B1 = 6.282509e3*cm^-1;
B3 = 1.9132288453;

C10 = 4.379437e12*cm^-1;
C12 = -1.314013e15*cm^-1;

% bVec = hpf(bVec);

%% Functions

% Short Range
USR = @(R) A1 + A2*(R/A0).^(-C);

% Intermediate Energy
xi = @(R) (R-Rm/A0)./(R+a*Rm/A0);
UIR = @ (Rvec) arrayfun(@(R) ...
    sum(bVec.*xi(R/A0).^(0:(length(bVec)-1))',1),Rvec);

% keyboard
% UIR = @ (R) UIR_func(R/A0,bVec,Rm/A0,a);



% Retardation Correction
fs = @(R) cos(R/lambda_bar) + ...
    (R/lambda_bar).*sin(R/lambda_bar);
fp = @(R) cos(R/lambda_bar) + ...
    (R/lambda_bar).*sin(R/lambda_bar)-...
    (R/lambda_bar).^2.*cos(R/lambda_bar);

% fs = @(R) 1;
% fp = @(R) 1;

% Long Range
Vp = @(R) D - fp(R).*C3p.*(R/A0).^(-3) ...
    - C6p*(R/A0).^(-6) - C8p*(R/A0).^(-8);
Vs = @(R) D - fs(R).*C3s.*(R/A0).^(-3) ...
    - C6s*(R/A0).^(-6) - C8s*(R/A0).^(-8);

% Adiabatic Energy
Vad_plus = @(R) -0.5*(Delta/3 - Vp(R) - Vs(R)) ...
    + 0.5*sqrt(8*(Delta/3)^2+(Vp(R)-Vs(R)-Delta/3).^2);
Vad_neg = @(R) -0.5*(Delta/3 - Vp(R) - Vs(R)) ...
    - 0.5*sqrt(8*(Delta/3)^2+(Vp(R)-Vs(R)-Delta/3).^2);



% Exchange
Vex = @(R) B1*(R/A0).^(B2).*exp(-B3*R/A0);

% Total long
ULR = @(R) Vad_neg(R)-Vex(R)-C10*(R/A0).^(-10)-C12*(R/A0).^(-12);

% Combined
U = @(R) (R<R1).*USR(R) + (R>R1).*(R<R2).*UIR(R) + (R>R2).*ULR(R);


%% 

Ra = linspace(2.55,R1/A0,1e3);
Rb = linspace(R1/A0,18.5,1e3);
Rc = linspace(18.5,50,1e3);




figure(22)
clf

plot(Ra,0.01*(USR(Ra*A0)-D),'r-','linewidth',2);
hold on
plot(Rb,0.01*(UIR(Rb*A0)-D),'b-','linewidth',2);
plot(Rc,0.01*(ULR(Rc*A0)-D),'g-','linewidth',2);
xlim([2.5 50])
plot(get(gca,'XLIm'),[1 1]*-0.01*2*Delta/3,'k--');
xlabel('separation (A)');
ylabel('energy cm^-1')

figure(20)
clf

subplot(211);
plot(Ra,0.01*(USR(Ra*A0)-D),'r-','linewidth',2);
hold on
plot(Rb,0.01*(UIR(Rb*A0)-D),'b-','linewidth',2);
plot(Rc,0.01*(ULR(Rc*A0)-D),'g-','linewidth',2);
xlim([15 25])
plot(get(gca,'XLIm'),[1 1]*-0.01*2*Delta/3,'k--');
xlabel('separation (A)');
ylabel('energy cm^-1')

subplot(212);
plot(Ra,0.01*(USR(Ra*A0)-D),'r-','linewidth',2);
hold on
plot(Rb,0.01*(UIR(Rb*A0)-D),'b-','linewidth',2);
plot(Rc,0.01*(ULR(Rc*A0)-D),'g-','linewidth',2);
xlim([2.5 3.5])
xlabel('separation (A)');
ylabel('energy cm^-1')