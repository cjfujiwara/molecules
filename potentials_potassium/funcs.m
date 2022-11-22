
% Intermediate Energy
xi = @(R) (R-Rm)./(R+b*Rm);
UIR = @ (Rvec) arrayfun(@(R) sum(aVec.*xi(R).^(0:(length(aVec)-1))',1),Rvec);

% Short Range
USR = @(R) A + (B/A0^Ns)*(R/A0).^(-Ns);

% Long Rangbe
Eexch = @(R) Aex*R.^gamma.*exp(-beta*R);

% Total Long Range Function
ULR = @(R) Uinf + ...
    -(C6/a0^6)*(R/a0).^(-6) + ...
    -(C8/a0^8)*(R/a0).^(-8) + ...
    -(C10/a0^10)*(R/a0).^(-10) + ...
    Eexch(R);

% Total Potential
U = @(R) (R<Rin).*USR(R) + (R>Rin).*(R<Rout).*UIR(R) + (R>Rout).*ULR(R);