function [engs,Ys,Xs]=solveRadial2(r,U,nMax)
%HYDROGEN_RADIAL Summary of this function goes here
Ufunc = @(rq) interp1(r,U,rq);

% Make differential equation for given eigenvalue and angular momentum
    function dydr = makeDiffEq(lambda,l)        
        dydr = @(r,y) [y(2);
            (Ufunc(r)+l*(l+1)/r^2-lambda)*y(1)];        
    end

% Calculate potential curve
lambda_vec=linspace(min(U)+.0000001,-.1,1000);
% doPlot = 0;
engs = zeros*ones(nMax,1);
errs = 100*ones(nMax,1);

ii = 1;
Xs={};
Ys={};
for n=1:length(lambda_vec)
    eng = lambda_vec(n);
    
    [n_nodes,err,Y,X]=findEig(eng);
    
    disp([num2str(n_nodes) ' ' num2str(err) ' ' num2str(eng)]);
    
    ii = n_nodes+1;
    try
        if abs(err)<abs(errs(ii))
           engs(ii) = eng;
           errs(ii) = err;
           Xs{ii} = X;
           Ys{ii} = Y;
        end
    end
    
    
end


 

function [n_nodes,err,Y,X]=findEig(lambda_0)
    % Find inner turning point and set it to be the matching point
    ind=find(U<lambda_0,1);
    r_L = r(ind)-1;

    % Find outer turning point and set it to be the matching point
    ind=find(flip(U<lambda_0),1);
    ind = length(U)-ind;    
    r_M = r(ind);

    % Outer calculation is larger than the turning point
    r_H = r_M+6;
    
    % Forward/Inner "L" Boundary Condition
    yf_0 = 0;
    zf_0 = 1;

    % Backward/Outer "H" Boundary Condition
    yb_inf = 0;
    zb_inf = -.001;

    % Make differential equation
    dydr = makeDiffEq(lambda_0,0);       

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

    err = yf(end)-yb(end);
end






   
end



