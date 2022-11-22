function output=solveRadial4(r,U,nMax,figstart)

%%
%% Settings

doPlot = 1;

%% Variables to Save
engs=zeros(nMax,1);
errs=zeros(nMax,1);
r_inner=zeros(nMax,1);
r_outer=zeros(nMax,1);
Xs={};
Ys={};

%% Initialize Strcuture
npt = struct;
% npt.r2U = @(rq) interp1(r,U,rq);

npt.r2U = @(rq) interp1qr(r,U,rq);
%% Calculate Curvature and HO energy

[~,imin]=min(U);
ii = [(imin-10):(imin+10)]';

pp=polyfit(r(ii),U(ii),2);
omega = 2*sqrt(pp(1));
dE = omega;
%% Find Ground State
% Initial n = 0
n = 1;
Ug = min(U) + dE*(0.5 + (n-1));

out = findState(r,U,Ug,n-1,0);

clear output
output(1) = out;

engs(n) = out.E;
errs(n) = out.err;
Xs{n} = out.X;
Ys{n} = out.Y;

str=[num2str(n-1) ' nodes : E=' num2str(out.E) ' e=' num2str(out.err)];
disp(str);

dEg = 2*(out.E - min(U));

if doPlot
    plotWavefunction(out,1+figstart);
end


%% Iterate on higher states

for n=2:nMax
    Ug = dEg + engs(n-1);
    
    if (nMax-n)<2
        doDetail = 1;
    else 
        doDetail = 0;
    end
    
    out =findState(r,U,Ug,n-1,doDetail);   

    output(n) = out;
    
    errs(n) = out.err;
    Xs{n} = out.X;
    Ys{n} = out.Y;
    engs(n) = out.E;    
    
    % Display output
    str=[num2str(n-1) ' nodes : E=' num2str(out.E) ' e=' num2str(out.err)];
    disp(str);

    % Find seperation to make new guess    
    if n<6    
        dEg = out.E-engs(n-1);
    else
        pp=polyfit((n-5):n,(engs((n-5):n))',2);
        Unew = polyval(pp,n+1);   
        dEg = Unew - engs(n);        
        if Unew < engs(n) || Unew > 0
            dEg = (out.E - engs(n-1))*0.1;
        end
    end 

    % Update Plots
    if doPlot
        plotWavefunction(out,figstart+n);       
    end

end

   
end

function plotWavefunction(wf,num)
    figure(99);
    clf
    set(gcf,'color','w','windowstyle','docked');
    co=get(gca,'colororder');
    clf
    
    yyaxis right
    p2=plot(wf.X,wf.Y(:,2),'--','color',co(2,:),'linewidth',1);     
%     hold on
    set(gca,'YColor',co(2,:)*.9);
    
    ylabel('derivative $y''(r)$','interpreter','latex');
    ym = max(abs(wf.Y(:,2)));
    ylim([-1 1]*ym);
    
    str=['$(n=' num2str(length(wf.r_node)) ',~E=' ...
        num2str(wf.E) ',\epsilon=' strrep(num2str(wf.err,'%.3E'),'E','\mathrm{E}') '$)'];    
    title(str,'interpreter','latex')
    
    set(gca,'box','on','linewidth',1,'fontname','times','fontsize',16);  
    
    yyaxis left
    p1=plot(wf.X,wf.Y(:,1),'-','color',co(1,:),'linewidth',2);
    
    xlim([min(wf.X) max(wf.X)]);
    ym = max(abs(wf.Y(:,1)));
    ylim([-1 1]*ym);
    hold on  
    plot(wf.r_node,zeros(length(wf.r_node),1),'o',...
        'markerfacecolor',co(1,:),'markersize',5,'markeredgecolor',co(1,:));    
    ylabel('wavefunction $y:=r R(r)$','interpreter','latex');    
    xlabel('separation ($\AA$)','interpreter','latex');

    plot([1 1]*wf.r_inner,get(gca,'YLim'),'k-');
    plot([1 1]*wf.r_outer,get(gca,'YLim'),'k-');
    
    set(gca,'YColor',co(1,:)*.9);
   drawnow;
end


function [r_L,r_M] = findTurningPoints(r,U,eng)
    r2U = @(rq) interp1(r,U,rq,'linear','extrap');

    r2U = @(rq) interp1qr(r,U,rq);

    
    ind=find(U<eng,1);
    r_L = r(ind);
    
    r_L = fzero(@(r) r2U(r)-eng,r_L);

    
    ind=find(flip(U<eng),1);
    ind = length(U)-ind;    
    r_M = r(ind);    
    r_M = fzero(@(r) r2U(r)-eng,r_M);  
end

function out = findState(r,U,Ug,n,doDetail)

out = struct;
%% SEtup input
npt = struct;
npt.E = Ug;
npt.r2U = @(rq) interp1(r,U,rq);
npt.r2U = @(rq) interp1qr(r,U,rq);

npt.n_node = n;
npt.U = U;

%% Adjust Bounds
% Find the classical turning points
[r_L_init,r_M_init] = findTurningPoints(r,U,Ug);

dR = r_M_init - r_L_init;

xL = 1/sqrt(-0.5*Ug);

r_H_init = r_L_init + dR*3;
r_L_init = r_L_init - 1;


if n >84
   r_H_init = r_L_init + dR*4;  
end

% if n==85
%    r_H_init = 250; 
% end

npt.r_L = r_L_init;
npt.r_M = r_M_init;
npt.r_H = r_H_init;

disp('initial integration');
[n_out,err,Y,X,x_nodes,Df,Db]=integrateSchrodinger(npt);

% Increase outer bounds
yd = abs(Y(:,2));
y = abs(Y(:,1));
increaseBounds = yd(end)/max(yd)>1e-4 || y(end)/max(y)>1e-3;
while increaseBounds 
    disp('increasing integration bounds');
    npt.r_H = npt.r_H + 2; 
    [n_out,err,Y,X,x_nodes,Df,Db]=integrateSchrodinger(npt);
    yd = abs(Y(:,2));
    increaseBounds = yd(end)/max(yd)>1e-4;
end


[~,npt.r_M] = findTurningPoints(r,U,Ug);
if n==85
   npt.r_H = 250;
end

%% Big Adjust


dN = n-n_out;
dm_sign = 0.5*(sign(Db)-sign(Df));
dm = Db - Df;


% if dN || abs(dm)>1
if dN || dm_sign
   warning(['Making adjustments ' num2str(n_out) ' nodes (' num2str(n) ...
        '); slope change sign : ' num2str(dm_sign)]);
    % Integrate schrodinger equation from one of the last few nodes. This
    % assumes that the integrated phase prior doesn't change much.
    % 
    % At the new rL y=0 and y'=c, where c sets an overal normalization. The
    % idea is to reduce the boundary substantially.
        
    n_out0 = n_out;
    err0=err;
    Y0=Y;
    X0=X;
    x_nodes0=x_nodes;
    Df0=Df;
    Db0=Db;
    
    % Spacing of last node
    dL_node = x_nodes0(end)-x_nodes0(end-1);
    
    % Position of last node
    r_node_last0 = x_nodes0(end);
    

    
    % Get temporary input
    nptSub = npt;
    nptSub.r_L = x_nodes(end-10);
    
    [n_nodes_s,err_s,Y_s,X_s,r_nodeVec_s,Df_s,Db_s]=integrateSchrodinger(nptSub);
    
    n_nodes_set = n_nodes_s + dN;    
    
    dE = min([1e-2 abs(npt.E)*0.05]);    
    if dN~=0
        warning(['Node Want : ' num2str(n) '; Node Have : ' num2str(n_out) '. Adjusting.']); 

        if dN > 0 
            sgn = +1;
        else
            sgn = -1;
        end
        
        % Iterate until N is correct.
        Niter = 0;
        
        while dN
            if nptSub.E+sgn*dE > 0
               dE=dE/2;            
            else
                nptSub.E = nptSub.E + sgn*dE;        
                fprintf([num2str(nptSub.E) ' ']);
                [n_nodes_s,err_s,Y_s,X_s,r_nodeVec_s,Df_s,Db_s]=integrateSchrodinger(nptSub);                
                
                % Get the last node
                r_node_last_new = r_nodeVec_s(end);     
                
                if (r_node_last_new > r_node_last0) && (abs(r_node_last_new-r_node_last0)>dL_node/2)
                   dN = dN - 1; 
                end
                
                if (r_node_last_new < r_node_last0) && (abs(r_node_last_new-r_node_last0)>dL_node/2)
                   dN = dN + 1; 
                end
            end      
            
            figure(20)
            clf
            plot(X_s,Y_s(:,1));
            title(Niter);
            
            Niter = Niter + 1;    

            if Niter>250
                
                keyboard
                return
            end  
        end
    end   
    
    % Recalculate matching point
    [~,r_M_init] = findTurningPoints(r,U,nptSub.E);
    nptSub.r_M = r_M_init;
    [n_nodes_s,err_s,Y_s,X_s,r_nodeVec_s,Df_s,Db_s]=integrateSchrodinger(nptSub);
    
    % Iterate until slopes have same sign    
    if sign(Db_s)~=sign(Df_s)
        warning('adjusting until slopes have same sign');
        if sign(Df_s)==1
            sgn = +1;   
        else
            sgn = -1;
        end
        Niter = 0;
        while sign(Df_s)~=sign(Db_s)
            if nptSub.E+sgn*dE > 0
                dE=dE/5;            
            else
                nptSub.E = nptSub.E + sgn*dE;        
                fprintf([num2str(nptSub.E) ' ']);
                [n_nodes_s,err_s,Y_s,X_s,r_nodeVec_s,Df_s,Db_s]=integrateSchrodinger(nptSub);
            end
            
            figure(20)
            clf
            plot(X_s,Y_s(:,1));
            Niter = Niter+1;
                        title(Niter);

        
            
            if Niter>100
                
                keyboard
                return
            end  
        end    
        

        warning('adjusting until slope magnitude changes');
        s0 = sign(Df_s-Db_s);
        s = s0;
        Niter=0;
        while s==s0
            if nptSub.E+sgn*dE > 0
                dE=dE/5;            
            else
                nptSub.E = nptSub.E + sgn*dE;        
                fprintf([num2str(nptSub.E) ' ']);
                [n_nodes_s,err_s,Y_s,X_s,r_nodeVec_s,Df_s,Db_s]=integrateSchrodinger(nptSub);
                s=sign(Df_s-Db_s);
            end  
            
                        
            figure(20)
            clf
            plot(X_s,Y_s(:,1));
                        title(Niter);

          Niter = Niter+1;
            
            if Niter>250
                
                keyboard
                return
            end  
        end     
    end                 
    npt.E = nptSub.E;    
end

%% Adjust matching point to turning point

[~,r_M_init] = findTurningPoints(r,U,npt.E);
npt.r_M = r_M_init;
[n_out,err,Y,X,x_nodes,Df,Db]=integrateSchrodinger(npt);   

if n_out~=n
   keyboard 
end


%% Reduce Bounds to improve integration time
% Find narrow bounds to improve integration time
P = Y(:,1).^2;
P = P / sum(P);

i1 = find(P/max(P)>1e-9,1);
i2 = length(P)-find(flip(P/max(P))>1e-9,1)+1;

r_L_new = X(i1);
r_H_new = X(i2);

npt.r_L = r_L_new;
npt.r_H = r_H_new;
npt.r_M = r_M_init;

if n_out ~= n 
    warning('oh no the nodal structure changed after changing bounds');
    keyboard
end
    
%%
    
% if n==85
%     npt.r_H = 250; 
% end
%% Optimize


Niter = 0;


Ug = npt.E;

eThresh = 0.04;
eThresh =  1e-2;

if abs(err) > eThresh
   fprintf('optimizing energy ...'); 
end


    dE = 0.01;
    dE = min([1e-2 abs(npt.E)*0.02]);    
while abs(err)> eThresh
    fprintf([num2str(npt.E) ',']);


    if Niter>100
        error('time out');
        return
    end
    
    Ug0 = npt.E;   
    n_out0 = n_out;
    err0 = err;
    Y0 = Y;
    X0 = X;
    nodes_0 = x_nodes;
    
    % Find error at displacement to find the slope
    
    if Ug0+dE>0
       dE = dE/10; 
    end
    
    npt.E = Ug0 + dE;
    [n_p,err_p,Yp,Xp,nodes_p]=integrateSchrodinger(npt);         
    Up = npt.E;
    dEdErr = dE/(err_p-err);     % Convert to slope
    
    % Modify
    Ug = Ug0 - 0.8*dEdErr*err0;   

    % Find updated error
     npt.E = Ug;
    [n_out,err,Y,X,x_nodes]=integrateSchrodinger(npt);  
    
    if n_out ~= n 
       warning('oh no the nodal structure changed');
       keyboard
    end
    
    Niter = Niter+1;
end

disp(' ');
Uout = Ug;

%% Final Update

if doDetail
%     [r_inner,r_outer] = findTurningPoints(r,U,npt.E);
%     npt.r_M = r_outer;
%     [n_out,err,Y,X,x_nodes]=integrateSchrodinger(npt);  
end


%% Find Turning Points

[r_inner,r_outer] = findTurningPoints(r,U,npt.E);


%%


out.E = Uout;
out.X = X;
out.Y = Y;
out.nodes = n_out;
out.r_inner = r_inner;
out.r_outer = r_outer;
out.r_node = x_nodes;
out.r_M = npt.r_M;
out.err = err;




end

function [n_nodes,err,Y,X,r_nodeVec,Df,Db]=integrateSchrodinger(npt)
    % Guess Energy
    E = npt.E;
    r2U = npt.r2U;
    
    r_L = npt.r_L;
    r_H = npt.r_H;
    r_M = npt.r_M;
    
    % Time independent schrodinger equation
    function dydr = makeDiffEq(lambda,l)        
        dydr = @(r,y) [y(2);
            (r2U(r)+l*(l+1)/r^2-lambda)*y(1)];        
    end


    
    % Forward/Inner "L" Boundary Condition
    yf_0 = 0;
    zf_0 = 1e-5;

    % Backward/Outer "H" Boundary Condition
    yb_inf = 0;
    zb_inf = -1e-120;

    % Make differential equation
    dydr = makeDiffEq(E,0);       

    options = odeset('refine',3);
    
    % Calculate integrationt
%     tic
    [rf,yf] = ode45(@(r,y) dydr(r,y),[r_L r_M],[yf_0;zf_0]);
%     toc
    [rb,yb] = ode45(@(r,y) dydr(r,y),[r_H r_M],[yb_inf;zb_inf]);   
%     toc
    
    if sum(sum(isnan(yb)))>0 || sum(sum(isnan(yf)))>0
        warning('NaN integration results')
        
    end
       
    % Match boundary conditions
    yf = yb(end,1)/yf(end,1)*yf;  

    % Find total wavefucntion and norm
    Y = [yf; flip(yb)]; X = [rf; flip(rb)];

    % Calculate number of nodes
%     p = Y(2:end-1,1);
%     rp = X(2:end-1,1);
    
%     p = Y(:,1);
%     rp = X(:,1);
    
%     p = yf(2:end-1,1);
%     rp = rf(2:end-1,1);

    p = yf(:,1);
    rp = rf(:,1);
    
    i_low = [abs(p)<1e-2*max(abs(p))];  
    p(i_low) = [];
    rp(i_low) = [];

    i_node = logical(abs(diff(sign(p)))*0.5);
    
    % Number of nodes
    n_nodes = sum(i_node);  
    
    % Find the position of the nodes by linearization
    ind_nodeVec = (1:length(i_node))';    
    ind_nodeVec = ind_nodeVec(i_node);
    r_nodeVec = zeros(length(ind_nodeVec),1);
    for kk=1:length(ind_nodeVec)
       r1 = rp(ind_nodeVec(kk));
       r2 = rp(ind_nodeVec(kk)+1);
       
       y1 = p(ind_nodeVec(kk));
       y2 = p(ind_nodeVec(kk)+1);
       
       m = (y2-y1)/(r2-r1);
       
       r_node = -y1/m + r1;
       
       r_nodeVec(kk) = r_node;
    end
    
    % Match boundary conditions
    y0 = max(Y(:,1));
    Y = Y/y0;
    yb = yb/y0;
    yf = yf/y0;

    % Calculate normalization
    N = sqrt(trapz(X,Y(:,1).^2));
    
    % Renormalize
    yf = yf/N;
    yb = yb/N; 
    Y = Y/N;
    
    if sum(sum(isnan(Y)))>1        
       warning('bad normalization or something');
       keyboard
    end

    % Calculate slope error
        
    Df=yf(end,2);
    Db=yb(end,2);
%     
    err = Df-Db;
%     keyboard
    err = (Df-Db)/max(abs([Df Db]));
    
%     err = yf(end,2)-yb(end,2);   

    
 
 
end




