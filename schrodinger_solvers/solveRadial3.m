function engs = solveRadial3(r,u)
l=0;
% Solve the radial Schrodinger equation with boundary conditions which has
% the form
% 
% d_rr(y) + [lambda - l(l+1)/r^2 - V(r)]y = 0 (eq. 1)
%
% Subject to the Dirchilet boundary conditions
%   y(0) = y(inf) = 0;
%
%% Transforming the Domain
% Modify equation (1) so that the domain is defined over the interval
% [-1,1] rather than [0,inf).
%
% x = (r-1)/(r+1)
% r = (1+x)/(1-x)

x2r = @(x) (1+x)./(1-x);
r2x = @(r) (r-1)./(r+1);


% r2u = @(rq) interp1(r,u,rq,'linear','extrap');
% 
% r2u = @(rq) max(Y(1),min(Y(end),interp1(X,Y,[0,1.5,4],'linear','extrap')))

%     function uout=r2u(rq)
%         i1 = rq<min(r);
%         i2 = rq>max(r);
%         
%         uout = interp1(r,u,rq,'linear','extrap');
%         uout(i1) = interp1(r,u,min(r),'linear','extrap');
%         uout(i2) = interp1(r,u,max(r),'linear','extrap');
%     end


r2u = @(rq) -2./rq;

x2u = @(x) r2u(x2r(x));

%
% So that at r=0 x=-1 and r=inf x=1.  With derivatives
%
% D[x,{r,1}] = 0.5*(x-1)^2
% D[x,{r,2}] = 0.5*(x-1)^3
%
% To change variables need to find d_xx(y) in terms of x
%
% dydr = dydx dxdr
% D[y,{r,2}] = D[y,{x,2}]*(D[x,{r,1}])^2 + D[y,{x,1}]*D[x,{r,2}]
% 
% d_rr(y) = 0.25*(x-1)^4*y''  + 0.5*(x-1)^3*y' 
%
% So that the new differential equation is
%
% 0.25*(x-1)^4*y''  + 0.5*(x-1)^3*y'  + [lambda - (1-x)^2/(1+x)^2*l(l+1) - V(x)]y = 0
%
% Or equivalently
%
% M1*D*D*y + M2*D*y + [-M3*l(l+1) - V[x]]y = -lambda*y
%% Grid Points
% To use Chebyshev, select an uneven grid which is provided by the cheb(N)
% function
%
N = 1000;

[D1,x] = cheb(N);

% Constructure second derivative matrix
D2 = D1*D1;

% % Enforce Boundary Conditions
% D1(1,:)   = zeros(1,N+1);
% D1(end,:) = zeros(1,N+1);
% 
% D2(1,:)   = zeros(1,N+1);
% D2(end,:) = zeros(1,N+1);


M1 = 0.25*diag((x-1).^4);       
M2 = 0.50*diag((x-1).^3);
M3 = diag((1-x).^2./(1+x).^2);


% M3(end,end) = 9999;

V = diag(x2u(x));
% V(1,1) = 0;
% 
% V(V==inf) = 9999;
% V(V==-inf) = -9999;


B = diag(-x.^2);
%% Assemble Matrix

H = M1*D2 + M2*D1  - V;

H(1,:) = zeros(1,N+1);
H(end,:) = zeros(1,N+1);
B(1,:) = zeros(1,N+1);
B(end,:) = zeros(1,N+1);
%% Diagonalize
[V,D] =  eig(H);

engs = sort(-diag(D));


end

% CHEB compute D = differentiation matrix, x = Chebyshev grid
function [D,x] = cheb(N)
    if N==0, D=0; x=1; return, end
    x = cos(pi*(0:N)/N)';
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X = repmat(x,1,N+1);
    dX = X-X';
    D = (c*(1./c)')./(dX+(eye(N+1))); % off-diagonal entries
    D = D - diag(sum(D')); % diagonal entries
end