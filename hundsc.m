%% Constants
Rvec = linspace(10,200,1000);
A0 = 10^(-10);

% Units
% energy = wavenumber * cL * h
% wavenumber = energy/(cL*h)

a0 = 5.29177210903e-11; % borh radius in m
e0 = 4.3597447222071e-18; % Hartree (energy) in J
cL = 299792458; % speed of light
h = 6.62607015e-34; % plancks' constant

amu = 1.66054e-27;

% m = 39.9637*amu;
m = 38.963706487*amu;

% Reduced mass
mu = m*m/(2*m);
%%
H0 = @(V1,V2,Delta) ...
    [V1-Delta/3 sqrt(2)/3*Delta;
    sqrt(2)/3*Delta V2];

H1 = @(V1,V2,V3,Delta) ...
    [V1 -Delta/3 Delta/3;
    -Delta/3 V2 Delta/3;
    Delta/3 Delta/3 V3];


H2 = @(V1) V1+Delta/3;

% Coefficients (in cm^-1)
C3p = 8.436*e0/(cL*h)*0.01;
C3s = 16.872*e0/(cL*h)*0.01;

C6p = 6272*e0/(cL*h)*0.01;
C6s = 9365*e0/(cL*h)*0.01;


C8p = 762300*e0/(cL*h)*0.01;
C8s = 1975000*e0/(cL*h)*0.01;

% Finestructure splitting
Delta = 57.7103143; % cm^-1
% Retardation Fcator

Lb = 767e-9/(2*pi);

fp = @(R) cos(R/Lb)+(R/Lb).*sin(R/Lb)-(R/Lb).^2.*cos(R/Lb);
fs = @(R) cos(R/Lb)+(R/Lb).*sin(R/Lb);

%% Create all data object

all_data = struct;

%% 0g+(3/2),0g+(1/2)
% 3Pig(rep) and 1sigma_g+(rep)


V1=@(R) +fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;
V2=@(R) +fs(R).*C3s./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;

engs = zeros(2,length(Rvec));

for kk=1:length(Rvec)
    v1 = V1(Rvec(kk)*A0);
    v2 = V2(Rvec(kk)*A0);
    H = H0(v1,v2,Delta);
    eng = eig(H);
    engs(:,kk) = eng;
end

all_data(1).Name = '0_g^+(1/2)';
all_data(1).Energy = engs(1,:);
all_data(1).LineStyle='--';

all_data(2).Name = '0_g^+(3/2)';
all_data(2).Energy = engs(2,:);
all_data(2).LineStyle='--';

%% 0u+(3/2),0u+(1/2)
% 3Piu(att) and 1sigma_u+(att)

V1=@(R) -fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;
V2=@(R) -fs(R).*C3s./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;

engs = zeros(2,length(Rvec));

for kk=1:length(Rvec)
    v1 = V1(Rvec(kk)*A0);
    v2 = V2(Rvec(kk)*A0);
    H = H0(v1,v2,Delta);
    eng = eig(H);
    engs(:,kk) = eng;
end

all_data(3).Name = '0_u^+(1/2)';
all_data(3).Energy = engs(1,:);
all_data(3).LineStyle='-';

all_data(4).Name = '0_u^+(3/2)';
all_data(4).Energy = engs(2,:);
all_data(4).LineStyle='-';


%% 0g-(3/2),0g-(1/2)
% 3Pig(rep) and 3sigma_g+(att)

V1=@(R) +fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;
V2=@(R) -fs(R).*C3s./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;

engs = zeros(2,length(Rvec));

dE2 = zeros(1,length(Rvec));

J = 0;


for kk=1:length(Rvec)
    v1 = V1(Rvec(kk)*A0);
    v2 = V2(Rvec(kk)*A0);
    H = H0(v1,v2,Delta);
    [b,eng] = eig(H);
    eng = diag(eng);
    engs(:,kk) = eng;
    
    
%     Erotp = hbar^2/(2*mu*(Rvec(kk)*A0)^2)*(J*(J+1)+2);
%     Erotp = Erotp/(cL*h*100);
% 
%     Erots = hbar^2/(2*mu*(Rvec(kk)*A0)^2)*(J*(J+1)+4);
%     Erots = Erots/(cL*h*100);
%     
%     dE2(kk)=b(1,2)*Erotp+b(2,2)*Erots;
end



all_data(5).Name = '0_g^-(1/2)';
all_data(5).Energy = engs(1,:);
all_data(5).LineStyle='-';

all_data(6).Name = '0_g^-(3/2)';
all_data(6).Energy = engs(2,:);
all_data(6).LineStyle='-';

%% 0u-(3/2),0u-(1/2)
% 3Piu(att) and 3sigma_u+(rep)

V1=@(R) -fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;
V2=@(R) +fs(R).*C3s./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;

engs = zeros(2,length(Rvec));

for kk=1:length(Rvec)
    v1 = V1(Rvec(kk)*A0);
    v2 = V2(Rvec(kk)*A0);
    H = H0(v1,v2,Delta);
    eng = eig(H);
    engs(:,kk) = eng;
end

all_data(7).Name = '0_u^-(1/2)';
all_data(7).Energy = engs(1,:);
all_data(7).LineStyle='--';

all_data(8).Name = '0_u^-(3/2)';
all_data(8).Energy = engs(2,:);
all_data(8).LineStyle='--';


%% 1g(3/2),1g(1/2),1g(3/2)
% 3Pig(rep) 1pig(att), 3sigmag+ (att)

V1=@(R) +fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;
V2=@(R) -fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;
V3=@(R) -fs(R).*C3s./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;

engs = zeros(3,length(Rvec));

for kk=1:length(Rvec)
    v1 = V1(Rvec(kk)*A0);
    v2 = V2(Rvec(kk)*A0);
    v3 = V3(Rvec(kk)*A0);

    H = H1(v1,v2,v3,Delta);
    eng = eig(H);
    engs(:,kk) = eng;
end

all_data(9).Name = '1_g(1/2)';
all_data(9).Energy = engs(1,:);
all_data(9).LineStyle='-';

all_data(10).Name = '1_g(3/2)';
all_data(10).Energy = engs(2,:);
all_data(10).LineStyle='-';

all_data(11).Name = '1_g(3/2)';
all_data(11).Energy = engs(3,:);
all_data(11).LineStyle='--';

%% 1u(3/2),1u(1/2),1u(3/2)
% 3Piu(att) 1piu(rep), 3sigmau+ (rep)

V1=@(R) -fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;
V2=@(R) +fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;
V3=@(R) +fs(R).*C3s./(R/a0).^3-C6p./(R/a0).^6-C8s./(R/a0).^8;

engs = zeros(3,length(Rvec));

for kk=1:length(Rvec)
    v1 = V1(Rvec(kk)*A0);
    v2 = V2(Rvec(kk)*A0);
    v3 = V3(Rvec(kk)*A0);

    H = H1(v1,v2,v3,Delta);
    eng = eig(H);
    engs(:,kk) = eng;
end

all_data(12).Name = '1_u(1/2)';
all_data(12).Energy = engs(1,:);
all_data(12).LineStyle='--';

all_data(13).Name = '1_u(3/2)';
all_data(13).Energy = engs(2,:);
all_data(13).LineStyle='-';

all_data(14).Name = '1_u(3/2)';
all_data(14).Energy = engs(3,:);
all_data(14).LineStyle='--';

%% 2g(3/2)
% 3Pig(rep)

V1=@(R) +fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;

all_data(15).Name = '2_g(3/2)';
all_data(15).Energy = V1(Rvec*A0)+Delta/3;
all_data(15).LineStyle='--';

%% 2u(3/2)
% 3Piu(att)

V1=@(R) -fp(R).*C3p./(R/a0).^3-C6p./(R/a0).^6-C8p./(R/a0).^8;

all_data(16).Name = '2_u(3/2)';
all_data(16).Energy = V1(Rvec*A0)+Delta/3;
all_data(16).LineStyle='--';

%% Groupings for labels

G1 = [3 5 9];
R1=19;

i1 = find(Rvec>R1,1);
e1 = zeros(length(G1),1);
for kk=1:length(G1)
    e1(kk) = all_data(G1(kk)).Energy(i1);
end
[~,inds]=sort(e1,'descend');
G1=G1(inds);


G2 = [16 10 12 4 7];
R2 = 19;

i2 = find(Rvec>R2,1);
e2 = zeros(length(G2),1);
for kk=1:length(G2)
    e2(kk) = all_data(G2(kk)).Energy(i1);
end
[~,inds]=sort(e2,'descend');
G2=G2(inds);


G3 = [1 6 13 11 15];
R3=19;

i3 = find(Rvec>R3,1);
e3 = zeros(length(G3),1);
for kk=1:length(G3)
    e3(kk) = all_data(G3(kk)).Energy(i1);
end
[~,inds]=sort(e3,'descend');
G3=G3(inds);


G4 = [2 14 8];
R4 = 19;

i4 = find(Rvec>R4,1);
e4 = zeros(length(G4),1);
for kk=1:length(G4)
    e4(kk) = all_data(G4(kk)).Energy(i1);
end
[~,inds]=sort(e4,'descend');
G4=G4(inds);


%%
hF=figure(10);
clf
hF.Color='w';
dco=get(gca,'colororder');
cmap = jet(16);

hF.Position = [10 100 800 600];

clf
plot([10 70],[1 1]*Delta/3,'k-','linewidth',1.5)
hold on
plot([10 70],[1 1]*-2*Delta/3,'k-','linewidth',1.5)

plot([1 1]*19,[-130 130],'k:','linewidth',1)

clear ps
for kk=1:length(all_data)
   ps(kk)=plot(Rvec,all_data(kk).Energy,all_data(kk).LineStyle,'Color',cmap(kk,:),...
       'linewidth',1);
   hold on
end

str1=[];


for kk=1:length(G1)
   ps(G1(kk)).Color=dco(4,:); 
   str1 = [str1 all_data(G1(kk)).Name newline];
end
str1(end)=[];

str2=[];
for kk=1:length(G2)
   ps(G2(kk)).Color=dco(1,:); 
  str2 = [str2 all_data(G2(kk)).Name newline];
end
str2(end)=[];

str3=[];
for kk=1:length(G3)
   ps(G3(kk)).Color=dco(3,:); 
     str3 = [str3 all_data(G3(kk)).Name newline];
end
str3(end)=[];

str4=[];
for kk=1:length(G4)
   ps(G4(kk)).Color=dco(2,:); 
    str4 = [str4 all_data(G4(kk)).Name newline];
end
str4(end)=[];


text(20,-150,str1,'edgecolor',dco(4,:),'backgroundcolor','w','fontsize',8,...
    'linewidth',2);
text(11,-55,str2,'edgecolor',dco(1,:),'backgroundcolor','w','fontsize',8,...
    'linewidth',2);
text(11,40,str3,'edgecolor',dco(3,:),'backgroundcolor','w','fontsize',8,...
    'linewidth',2);
text(20,130,str4,'edgecolor',dco(2,:),'backgroundcolor','w','fontsize',8,...
    'linewidth',2);

s1 = '${^2S_{1/2}}+{^2 P_{3/2}}$';
s2 = '${^2S_{1/2}}+{^2 P_{1/2}}$';
s3 = '$\mathrm{K}_2$ Long Range';

text(39.5,60,s1,'horizontalalignment','right','interpreter','latex',...
    'fontsize',18);
text(39.5,-70,s2,'horizontalalignment','right','interpreter','latex',...
    'fontsize',18);

text(.98,.98,s3,'interpreter','latex','units','normalized','fontsize',24,...
    'horizontalalignment','right','verticalalignment','cap');

xlim([10 40]);

ylim([-200 200]);

legend(ps,{all_data.Name},'location','eastoutside')

set(gca,'xgrid','on','ygrid','on','yminorgrid','on','fontsize',12);

ylabel('energy (cm^-1)');
xlabel(['separation (' char(197) ')'])

%% Zoom in on Long Range Portion
hF=figure(11);
clf
hF.Color='w';
dco=get(gca,'colororder');
cmap = jet(16);

hF.Position = [10 100 800 600];

clf
plot([10 70],[1 1]*0,'k-','linewidth',1.5)
hold on


clear ps
for kk=1:length(all_data)
   ps(kk)=plot(Rvec,cL*100*1e-9*(all_data(kk).Energy-Delta/3),all_data(kk).LineStyle,'Color',cmap(kk,:),...
       'linewidth',2);
   hold on
end

str1=[];


for kk=1:length(G1)
   ps(G1(kk)).Color=dco(4,:); 
   str1 = [str1 all_data(G1(kk)).Name newline];
end
str1(end)=[];

str2=[];
for kk=1:length(G2)
   ps(G2(kk)).Color=dco(1,:); 
  str2 = [str2 all_data(G2(kk)).Name newline];
end
str2(end)=[];

str3=[];
for kk=1:length(G3)
   ps(G3(kk)).Color=dco(3,:); 
     str3 = [str3 all_data(G3(kk)).Name newline];
end
str3(end)=[];

str4=[];
for kk=1:length(G4)
   ps(G4(kk)).Color=dco(2,:); 
    str4 = [str4 all_data(G4(kk)).Name newline];
end
str4(end)=[];

xlim([22 50]);

ylim([-9 3]*cL*100/1e9);

% Strings for the specific states
ss={'$1_u(3/2)$','$0_g^-(3/2)$','$2_u(3/2)$','$1_g(3/2)$','$0_u^+(3/2)$'};

x=[35 26 28.5 36 38.5];
y=[-30 -210 -250 -200 -250];

for kk=1:length(ss)
    text(x(kk),y(kk),ss{kk},'interpreter','latex','fontsize',16);
end



s1 = '${^2S_{1/2}}+{^2 P_{3/2}}$';
text(40,12,s1,'horizontalalignment','right','interpreter','latex',...
    'fontsize',24);

set(gca,'xgrid','on','ygrid','on','yminorgrid','on','fontsize',12);

ylabel('energy (GHz)');
xlabel(['separation (' char(197) ')'])

%% Diagnolize long range potentials

me = 9.1093837015e-31; % electron mass [kg]
hbar = 1.054571817e-34; % reduced bplanck's constant [Js]
h = hbar*(2*pi); % plcnak's constant
Eh = hbar^2/(me*a0^2);  % Hartree energy [J]
c = 299792458; % speed of light [m/s]

%% 0g-(3/2)


r = Rvec;


% 0g-(3/2)
i1 = 6;
u1_cm = (all_data(i1).Energy-Delta/3);
u1_blah = (100*u1_cm*cL*h)/(hbar^2/(2*mu*A0^2));


u_test_cm = -6.4933+0.24870;
u_test_blah = (u_test_cm*100*cL*h)/(hbar^2/(2*mu*A0^2));
% solveRadial(r,u1_blah,u_test_blah)

nMax=31;

[engs_out,Ys,Rs]=solveRadial2(r,u1_blah,nMax);

%% 1g(3/2)

% 
% r = Rvec;
% 
% 
% % 1g(3/2)
% i1 = 10;
% u1_cm = (all_data(i1).Energy-Delta/3);
% u1_blah = (100*u1_cm*cL*h)/(hbar^2/(2*mu*A0^2));
% 
% 
% u_test_cm = -6.4933+0.24870;
% u_test_blah = (u_test_cm*100*cL*h)/(hbar^2/(2*mu*A0^2));
% % solveRadial(r,u1_blah,u_test_blah)
% 
% nMax=31;
% 
% [engs_out,Ys,Rs]=solveRadial2(r,u1_blah,nMax);

%%
engs_detune_GHz = 1e-9*engs_out * (hbar^2/(2*mu*A0^2)) / (h);
engs_cm = (engs_out-min(u1_blah)) * (hbar^2/(2*mu*A0^2)) / (h*100*cL);

summary = [engs_detune_GHz engs_cm];

disp(' ');
disp(summary);

%% Zoom in on Long Range Portion
hF=figure(12);
clf
hF.Color='w';
dco=get(gca,'colororder');
cmap = jet(16);

hF.Position = [10 100 800 600];
clf

xlim([22 50]);
i1 = 6;
u1_cm = (all_data(i1).Energy-Delta/3);
u1_GHz = 100*u1_cm*cL/1e9;
plot(r,u1_GHz,'linewidth',2,'color','r');
hold on
for kk=1:31
    u_n_GHz = engs_detune_GHz(kk);
    i1 = find(u_n_GHz>u1_GHz,1);
    i2 = find(u_n_GHz>flip(u1_GHz),1);
    i2 = length(u1_GHz)-i2;    
    plot([r(i1)-.2 r(i2)+.5],[1 1]*u_n_GHz,'k-');
end

ylim([- 200 0]);
xlim([20 50]);

set(gca,'xgrid','on','ygrid','on','yminorgrid','on','fontsize',12);
ylabel('energy (GHz)');
xlabel(['separation (' char(197) ')'])
title('39$\mathrm{K}_2$ $0_g^-$ vibrational modes','interpreter','latex');

