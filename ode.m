function dYdt = ode(t,Y,M)
Uold = Y(1:M.N);Vold = Y(M.N+1:2*M.N);Qold = Y(2*M.N+1:3*M.N);

% Calculate the effective pressure
P = NormalStress(t,M);

% Friction
Tfriction = M.a .* P .* asinh( Vold/M.V0/2 .* exp( Qold ./ M.a ));

% Floating ice has no friction
if numel(M.floating)>0, Tfriction(M.floating) = 0; end

% Tides
Tfront = TideStress(t,M);

% Shear
if M.N > 1
    Txy_y = M.G * M.D2 * Uold;
else
    Txy_y = 0;
end

% Upstream loading
Txx = M.G*2*(1-M.poi)/(1-2*M.poi) * (M.Vup*t-Uold) / M.Lup;

% Field Updates
dU =  Vold; % Only used for elastic case
dV =  (Txy_y - (Tfriction-M.Tau0)/M.H + Txx/M.Lup + Tfront/M.Ldown)/M.rho;

% Slip Law
dQ = -abs(Vold)./M.Dc .* (Qold - M.f0 + M.b.*log(abs(Vold)./M.V0));

% Ageing law
% dQ =  M.b.*M.V0./M.Dc.* ( exp( (M.f0-Qold)./M.b ) - abs(Vold)./M.V0 );

% Floating ice has no state evolution
if numel(M.floating)>0, dQ(M.floating) = 0; end

dYdt = [dU;dV;dQ];