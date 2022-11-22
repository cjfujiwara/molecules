function hydrogen_radial
%HYDROGEN_RADIAL Summary of this function goes here

% Hydrodgen potential
V = @(r) -2./r;

% Make differential equation for given eigenvalue and angular momentum
    function dydr = makeDiffEq(lambda,l)        
        dydr = @(r,y) [y(2);
            (V(r)+l*(l+1)/r^2-lambda)*y(1)];        
    end


% Matching points
r_H = 30;
r_L = 1e-9;

% Forward/Inner "L" Boundary Condition
yf_0 = 0;
zf_0 = 1;

% Backward/Outer "H" Boundary Condition
yb_inf = 0;
zb_inf = -.001;

% Calculate potential curve
xx=linspace(0,100,1000);
Ur=-2./xx;

lambda_vec=linspace(-1,-.0278,300);

lambda_vec=-1./(linspace(1,6,300)).^2;

for kk=1:length(lambda_vec)
    % Energy Guess
%     lambda_0=-.111;
    lambda_0 = lambda_vec(kk);
    
    % Find outer turning point and set it to be the matching point
    ind=find(Ur>lambda_0,1);
    r_M = xx(ind);
    r_H = 3*r_M;
    
    % Make differential equation
    dydr = makeDiffEq(lambda_0,0);

    % Calculate integration
    [rf,yf] = ode45(@(r,y) dydr(r,y),[r_L r_M],[yf_0;zf_0]);
    [rb,yb] = ode45(@(r,y) dydr(r,y),[r_H r_M],[yb_inf;zb_inf]);

    % Match boundary conditions
    yf = yb(end,1)/yf(end,1)*yf;

    % Find total wavefucntion and norm
    Y = [yf; flip(yb)]; X = [rf; flip(rb)];
    N = sqrt(trapz(X,Y(:,1).^2));

    % Renormalize
    yf = yf/N;
    yb = yb/N;

    hF=figure(10);
    clf
    hF.Color='w';
    
    subplot(211);
    title('Coulomb Radial Eigenfunctions');

    plot(rf,yf(:,1),'r-');
    hold on
    plot(rb,yb(:,1),'b-');
    xlabel('\rho');
    ylabel('y(\rho)');
    xlim([0 r_H]);
    yyaxis right

    plot(xx,Ur,'k:');
    hold on
    plot([0 r_H],[1 1]*lambda_0,'k-');
    xlim([0 r_H]);


    ylim([-1.2 0]);
    set(gca,'YColor','k');
    ylabel('energy');
    xlim([0 r_H]);
    
    str = ['$\lambda = -1/(' num2str(1./sqrt(-lambda_0)) ')^2$'];
    text(.98,.02,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'horizontalalignment','right',...
        'verticalalignment','bottom');

    subplot(212);
    plot(rf,yf(:,2),'r-');
    hold on
    plot(rb,yb(:,2),'b-');
    xlabel('\rho');
    ylabel('y''(\rho)');
    xlim([0 r_H]);
    

%     F(kk) = getframe(gcf) ;
    drawnow;
end



%     
%   % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 10;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);


end

