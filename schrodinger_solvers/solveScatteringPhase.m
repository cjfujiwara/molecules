function [delta] = solveScatteringPhase(r,U,Escat,l)
q = sqrt(Escat);

if max(r)<2*2*pi/q
    warning('not many periods for scattering length to be in');
end


%% Integrate schrodinger equation
r1=findInnerTurningPoint(r,U,Escat);

npt = struct;
npt.r2U = @(rq) interp1qr(r,U,rq);
npt.r_L = r1-.2;
npt.r_H = max(r);
npt.E = Escat;



[X,Y]=integrateSchrodinger(npt);

i = find(X>200,1);

Xl = X(i:end);
Yl = Y(i:end,1);

% Find location of first node
s = sign(Yl);
ds = diff(s);
n1 = ds*.5;
i1 = find(n1==1,1);


% Location of a node with - slope
x1 = Xl(i1);
delta0 = mod(-x1*q,2*pi);
A0 = max(Yl);

%% Construct Fit

foo = @(A,delta,rho) A*sin(q*rho - pi*l/2+ delta);

myfit = fittype(@(A,delta,rho) foo(A,delta,rho),'independent',{'rho'},...
    'coefficients',{'A','delta'});
fitopt = fitoptions(myfit);
fitopt.StartPoint = [A0 delta0];

fout=fit(Xl,Yl,myfit,fitopt);

delta = fout.delta;
end

% Calculate the classical turning points
function [r_L] = findInnerTurningPoint(r,U,eng)

    % Interpolant function
    r2U = @(rq) interp1qr(r,U,rq);

    % Find the crossing point
    ind=find(U<eng,1);
    
    % Numerically find the actual zero crossing
    r_L = r(ind);    
    r_L = fzero(@(r) r2U(r)-eng,r_L);

end

function [X,Y]=integrateSchrodinger(npt)
    % Guess Energy
    E = npt.E;
    r2U = npt.r2U;
    
    
    r_L = npt.r_L;
    r_H = npt.r_H;
    
    % Time independent schrodinger equation
    function dydr = makeDiffEq(lambda,l)        
        dydr = @(r,y) [y(2);
            (r2U(r)+l*(l+1)/r^2-lambda)*y(1)];        
    end
    
    % Forward/Inner "L" Boundary Condition
    yf_0 = 0;
    zf_0 = 1e-5;

    % Make differential equation
    dydr = makeDiffEq(E,0);       

    options = odeset('refine',5);
    
    % Calculate integrationt
    [X,Y] = ode45(@(r,y) dydr(r,y),[r_L r_H],[yf_0;zf_0],options);
    
    if sum(sum(isnan(Y)))>0
        warning('NaN integration results')        
    end      
 
end
