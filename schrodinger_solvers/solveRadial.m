function solveRadial(r,U,lambda_0)
%HYDROGEN_RADIAL Summary of this function goes here
Ufunc = @(rq) interp1(r,U,rq);

doPlot = 1;

% Make differential equation for given eigenvalue and angular momentum
    function dydr = makeDiffEq(lambda,l)        
        dydr = @(r,y) [y(2);
            (Ufunc(r)+l*(l+1)/r^2-lambda)*y(1)];        
    end

% Forward/Inner "L" Boundary Condition
yf_0 = 0;
zf_0 = 1;

% Backward/Outer "H" Boundary Condition
yb_inf = 0;
zb_inf = -.001;




% Find inner turning point and set it to be the matching point
ind=find(U<lambda_0,1);
r_L = r(ind)-1;

% Find outer turning point and set it to be the matching point
ind=find(flip(U<lambda_0),1);
ind = length(U)-ind;    
r_M = r(ind);

% Outer calculation is larger than the turning point
r_H = r_M+3;

% Make differential equation
dydr = makeDiffEq(lambda_0,0);   

disp(r_L)
disp(r_H);
disp(r_M);

% Calculate integration
[rf,yf] = ode45(@(r,y) dydr(r,y),[r_L r_M],[yf_0;zf_0]);
[rb,yb] = ode45(@(r,y) dydr(r,y),[r_H r_M],[yb_inf;zb_inf]);

% Match boundary conditions
yf = yb(end,1)/yf(end,1)*yf;

% Find total wavefucntion and norm
Y = [yf; flip(yb)]; X = [rf; flip(rb)];

n_nodes = 0.5*sum(abs(diff(sign(Y(2:end-1,1)))));

N = sqrt(trapz(X,Y(:,1).^2));

% Renormalize
yf = yf/N;
yb = yb/N;

if doPlot
    hF=figure(14);
    clf
    hF.Color='w';

    subplot(211);
    title('Coulomb Radial Eigenfunctions');

    plot(rf,yf(:,1),'r-');
    hold on
    plot(rb,yb(:,1),'b-');
    xlabel('\rho (a_0)');
    ylabel('y(\rho)');
    xlim([r_L r_H]);
    yyaxis right

    plot(r,U,'k:');
    ylim([-1 1]+lambda_0);
    hold on
    plot([r_L r_H],[1 1]*lambda_0,'k-');
    xlim([r_L r_H]);

    plot([1 1]*r_M,get(gca,'YLim'),'k-');

    set(gca,'YColor','k');
    ylabel('energy');
    xlim([r_L r_H]);

    str = ['$\lambda = -1/(' num2str(1./sqrt(-lambda_0)) ')^2$'];
    text(.98,.02,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'horizontalalignment','right',...
        'verticalalignment','bottom');

    str = [num2str(n_nodes) ' nodes'];
    text(.02,.05,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'horizontalalignment','left',...
        'verticalalignment','bottom');

    subplot(212);
    plot(rf,yf(:,2),'r-');
    hold on
    plot(rb,yb(:,2),'b-');
    xlabel('\rho (a_0)');
    ylabel('y''(\rho)');
    xlim([r_L r_H]);
end
    
    
   
end




   



