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
RR = linspace(1,300,2e4)';

% Potential for X^1 sigma^+_g (ground state)
K_X_1_sigma;
funcs;
U_X_1_sigma = U(RR*A0);

% Potential for a^3 sigma^+_u (first excited)
K_a_3_sigma;
funcs;
U_a_3_sigma = U(RR*A0);


%% Plot Potential Curves

hF=figure(10);
hF.Position=[100 500 600 400];
clf
set(gcf,'color','w');

% Zoom out potential energies
ax1=axes;
p1=plot(RR,U_X_1_sigma*0.01/10^3,'linewidth',2);
hold on
p2=plot(RR,U_a_3_sigma*0.01/10^3,'linewidth',2);

% X and Y limits and other labels
xlim([2.5 14.1]);
ylim([-4.5 .8]);
ylabel('V (10^3 cm^-1)');
xlabel(['separation (' char(197) ')']);
set(gca,'xgrid','off','ygrid','off','fontsize',14,...
    'box','on','linewidth',1,'units','normalized',...
    'fontname','times');

% Match right boundary to be in THz
yL = get(gca,'YLim');
yL_m=yL*1e5; % convert to inverse meters
yL_Hz = yL_m*cL; % Limits in Hz
yyaxis right
ylim(yL_Hz*1e-12);
ylabel('energy (THz)');
set(gca,'YColor','k');
yyaxis left

% Create for plots
legStr = {'$X~^{1}\Sigma^+_g$','$a~^{3}\Sigma^+_u$'};
legend([p1 p2],legStr,'interpreter','latex','fontsize',14,...
    'location','northeast','orientation','horizontal');

% Secondary axis for zoom in
ax2 = axes('units','normalized');
ax2.Position(3)=.25;
ax2.Position(4)= .2;
ax2.Position(1) = ax1.Position(1)-ax2.Position(3)+ax1.Position(3);
ax2.Position(2) = ax1.Position(2);

% Plot on the zoom in
plot(RR,U_X_1_sigma*.01,'linewidth',2);
hold on
plot(RR,U_a_3_sigma*.01,'--','linewidth',2);

% Label the zoom in
xlim([10.5 14.5]);
ylim([-17 0]);
set(gca,'ytick',[-15 -10 -5 0],'xaxislocation','top',...
    'xtick',[11 12 13 14]);
ylabel('V $(\mathrm{cm}^{-1})$','interpreter','latex','fontsize',10);
xlabel(['separation ($\AA$)'],'interpreter','latex');


%% Diagonalize X state
U_xsigma_blah = (U_X_1_sigma*cL*h)/(hbar^2/(2*mu*A0^2));
if m==39
    nMax = 86;
else 
    nMax = 87;
end

Xsigma=solveRadial4(RR,U_xsigma_blah,nMax,200);

% Convert units to useful ones
eng_J = [Xsigma.E] * (hbar^2/(2*mu*A0^2));
eng_cm = eng_J/(cL*h);
eng_cm = num2cell(eng_cm);
[Xsigma(:).E_cm] = deal(eng_cm{:});

%% Diagonalize a State
if m==39
    nMax = 27;
else 
    nMax = 27;
end

U_asigma_blah = (U_a_3_sigma*cL*h)/(hbar^2/(2*mu*A0^2));

figstart=100;
asigma=solveRadial4(RR,U_asigma_blah,nMax,figstart);

% Convert units to useful ones
eng_J = [asigma.E] * (hbar^2/(2*mu*A0^2));
eng_cm = eng_J/(cL*h);
eng_cm = num2cell(eng_cm);
[asigma(:).E_cm] = deal(eng_cm{:});

%% Plot the Energies vs eigenvalue
hF=figure(2);
hF.Position=[100 100 600 400];

co = get(gca,'colororder');
clf
set(gcf,'color','w');

% X1 energies
X_Ub = min(U_X_1_sigma);
X_Eb13 = (-[Xsigma.E_cm]).^(1/3);
X_nu = [Xsigma.nodes];
X_pp=polyfit(X_nu(end-6:end),X_Eb13(end-6:end),1);

X_slope_str = [num2str(X_pp(1)) ' $\mathrm{cm}^{1/3}$'];

% a3 energies
a_Ub = min(U_a_3_sigma);
a_Eb13 = (-[asigma.E_cm]).^(1/3);
a_nu = [asigma.nodes];
a_pp=polyfit(a_nu(end-6:end),a_Eb13(end-6:end),1);

a_slope_str = [num2str(a_pp(1)) ' $\mathrm{cm}^{1/3}$'];

tt=linspace(0,90,10);

plot(X_nu,X_Eb13,'o','markerfacecolor',co(1,:),'markersize',6,...
    'linewidth',2,'markeredgecolor',co(1,:)*.8);
% xlabel('$X^1\Sigma_g^+$ mode number','interpreter','latex');
ylabel('binding energy $E_b^{1/3}~(\mathrm{cm}^{-1/3})$','interpreter','latex');
xlim([-6 0]+max(X_nu)+1);
ylim([0 10]);
set(gca,'box','off','linewidth',1,'fontname','times','fontsize',14,...
    'XColor',.9*co(1,:));
hold on
pxf=plot(tt,X_pp(1)*tt+X_pp(2),'-','color',co(1,:)*.8,'linewidth',2);


hAx(1)=gca;
hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
hold(hAx(2),'on')
pa=plot(hAx(2),a_nu,a_Eb13,'o','markerfacecolor',co(2,:),'markersize',6,...
    'linewidth',2,'markeredgecolor',co(2,:)*.8);
set(gca,'box','off','linewidth',1,'fontname','times','fontsize',14,...
    'xgrid','on','ygrid','on','XColor',.9*co(2,:));
xlim([-6 0]+max(a_nu)+1);
ylim([0 10]);
% xlabel('$a^3\Sigma_u^+$ mode number','interpreter','latex');
hold on
paf=plot(tt,a_pp(1)*tt+a_pp(2),'-','color',co(2,:)*.8,'linewidth',2);


text(hAx(1),.5,.02,'$X^1\Sigma_g^+$ mode number','interpreter','latex',...
    'units','normalized','horizontalalignment','center',...
    'verticalalignment','bottom','fontsize',20,'color',co(1,:)*.8);
text(hAx(2),.5,.98,'$a^3\Sigma_u^+$ mode number','interpreter','latex',...
    'units','normalized','horizontalalignment','center',...
    'verticalalignment','top','fontsize',20,'color',co(2,:)*.9);

text(hAx(2),.98,.98,mass_str,'interpreter','latex',...
    'units','normalized','horizontalalignment','right',...
    'verticalalignment','top','fontsize',24,'color','k');

legend([paf pxf],{a_slope_str,X_slope_str},'interpreter','latex',...
    'location','east');
%% Plot Energies All

hF=figure(12);
hF.Position=[100 100 900 500];
clf
set(gcf,'color','w');
co=get(gca,'colororder');
clf;


% Xsigma
subplot(121)

% plot potential
plot(RR,U_X_1_sigma*0.01/10^3,'linewidth',2,'color',co(1,:));

hold on

% Plot eigenvalues
for kk=1:length(Xsigma)    
    E = Xsigma(kk).E_cm;
    r_1 = Xsigma(kk).r_inner;
    r_2 = Xsigma(kk).r_outer;    
    plot([r_1 r_2],[1 1]*E*0.01/10^3,'-','color',co(1,:),'linewidth',1); 
end

% Refining labels
xlim([2.5 11]);
ylim([-4.5 .05]);
ylabel('V (10^3 cm^-1)');
xlabel(['separation (' char(197) ')']);
str = {'$X~^{1}\Sigma^+_g$'};
legend(str,'interpreter','latex','fontsize',16,'location','southeast');
set(gca,'xgrid','on','ygrid','on','fontsize',12,...
    'box','on','linewidth',1,'fontname','times');

text(.95,.2,mass_str,'interpreter','latex',...
    'units','normalized','horizontalalignment','right',...
    'verticalalignment','top','fontsize',24,'color','k');

% asigma
subplot(122)
plot(RR,U_a_3_sigma*0.01/10^3,'linewidth',2,'color',co(2,:));
hold on

% Plot eigenvalues
for kk=1:length(asigma)    
    E = asigma(kk).E_cm;
    r_1 = asigma(kk).r_inner;
    r_2 = asigma(kk).r_outer;    
    plot([r_1 r_2],[1 1]*E*0.01/10^3,'-','color',co(2,:),'linewidth',1); 
end

% Refining labels
xlim([2.5 11]);
ylim([-.26 .005]);
ylabel('V (10^3 cm^-1)');
xlabel(['separation (' char(197) ')']);
str={'$a~^{3}\Sigma^+_u$'};
set(gca,'xgrid','on','ygrid','on','fontsize',12,...
    'box','on','linewidth',1,'fontname','times');
legend(str,'interpreter','latex','fontsize',16,'location','southeast');



text(.95,.2,mass_str,'interpreter','latex',...
    'units','normalized','horizontalalignment','right',...
    'verticalalignment','top','fontsize',24,'color','k');
%% Zoom in on final states

str_vdw = '$R_\matrhm{vdW} = \frac{1}{2}\left(\frac{2\mu C_6}{\hbar^2}\right)^4$';
str_evdw = '$E_\mathrm{vdW} = \frac{\hbar^2}{2\mu}\frac{1}{R_\mathrm{vDW}^2}$';


Rvdw = 0.5*(2*mu*C6*cL*h /hbar^2)^(1/4); % van der waals length [m]
Evdw = (hbar^2/(2*mu))*(1/Rvdw^2)/h*1e-9; % van der waals length [Ghz]

hF=figure(13);
co=get(gca,'colororder');
hF.Position=[100 100 600 600];
clf
set(gcf,'color','w');

% Plot potentials (versus bohr)
plot([1 1]*Rvdw/a0,[-3 2],'k--');
% plot([0 1]*Rvdw/a0,-[1 1]*Evdw,'k-');

hold on

p1=plot(RR*A0/a0,U_X_1_sigma*cL*1e-9,'-','linewidth',3,'color',co(1,:));
hold on
p2=plot(RR*A0/a0,U_a_3_sigma*cL*1e-9,'-','linewidth',3,'color',co(2,:));
legStr = {'$X{^{1}\Sigma}^+_g$','$a{^{3}\Sigma}^+_u$'};

hold on


% Plot eigenvalues

% L_extra = zeros(
for kk=(length(Xsigma)-1):length(Xsigma)    
    
    E = Xsigma(kk).E_cm*cL*1e-9;
    r_1 = Xsigma(kk).r_inner*A0/a0;
    r_2 = Xsigma(kk).r_outer*A0/a0;
    
    Xe = Xsigma(kk).X(end)*A0/a0;
    plot([r_1 r_2],[1 1]*E,'-','color',co(1,:),'linewidth',2); 
    plot([r_2 Xe+5],[1 1]*E,'--','color',co(1,:),'linewidth',1); 
    
    str = ['$\nu=' num2str(Xsigma(kk).nodes) '$' newline ...
        '$' num2str(round(1e3*E,1)) '~\mathrm{MHz},~' num2str(round(r_2,1)) 'a_0$'];

    text(195,E+.02,str,'verticalalignment',...
        'bottom','fontsize',20,'fontname','times',...
        'color',co(1,:),'interpreter','latex',...
        'horizontalalignment','right'); 
    
    A = .2;
    Y = Xsigma(kk).Y(:,1);
    Y = (A/max(Y))*Y;
    X = Xsigma(kk).X*A0/a0;
    
    plot(X,Y+E,'-','color',co(1,:),'linewidth',1);
    
end

% Plot eigenvalues
for kk=(length(asigma)-1):length(asigma)    
    E = asigma(kk).E_cm*cL*1e-9;
    r_1 = asigma(kk).r_inner*A0/a0;
    r_2 = asigma(kk).r_outer*A0/a0;        
    Xe = asigma(kk).X(end)*A0/a0;

    plot([r_1 r_2],[1 1]*E,'-','color',co(2,:),'linewidth',2); 
    plot([r_2 Xe+5],[1 1]*E,'--','color',co(2,:),'linewidth',1); 

    
    str = ['$\nu=' num2str(kk) '$' newline ...
        '$' num2str(round(1e3*E,1)) '~\mathrm{MHz},~' num2str(round(r_2,1)) 'a_0$'];

    text(195,E+.02,str,'verticalalignment',...
        'bottom','fontsize',20,'fontname','times',...
        'color',co(2,:),'interpreter','latex',...
        'horizontalalignment','right'); 
    
    
    A = .2;
    Y = asigma(kk).Y(:,1);
    Y = (A/max(Y))*Y;
    X = asigma(kk).X*A0/a0;
    
    plot(X,Y+E,'-','color',co(2,:),'linewidth',1);
    
end


xlim([2.5 200]);
ylim([-2.3 .5]);

ylabel('energy (GHz)','interpreter','latex');
% xlabel(['separation (' char(197) ')']);
xlabel(['separation ($a_0$)'],'interpreter','latex');

str = {'$X~^{1}\Sigma^+_g$'};


legend([p1 p2],legStr,'interpreter','latex','fontsize',18,'location','southeast',...
    'orientation','horizontal');



set(gca,'xgrid','on','ygrid','on','fontsize',16,...
    'box','on','linewidth',1,'fontname','times');

text(.05,.98,mass_str,'interpreter','latex',...
    'units','normalized','horizontalalignment','left',...
    'verticalalignment','top','fontsize',24,'color','k');
