function [Dx,A,A2] = createSBPOps_mc(Nx,dx,cIn,order)
% creates narrow-stencil variable coefficient matrices in 2D
% matrices are compatible

c = sparse(diag(cIn(:))); % coefficient matrix


% 1D matrix factors
xMats = sbp1d_mc(Nx,dx,order);


Dx = xMats.D1; % 1st derivative



% create 2nd derivative with variable coefficients
% (This formula comes from equation 8 in
% http://link.springer.com/article/10.1007/s10915-011-9525-z)
if order==2
  C2x = speye(Nx,Nx);C2x(1,1)=0;C2x(Nx,Nx)=0;
  Rxmu = xMats.D2'*C2x*c*xMats.D2.*(1/4/dx);
elseif order==4
  cc = cIn(:);
  B3 = 0.5*sparse(1:Nx,1:Nx,cc,Nx,Nx)...
    + 0.5*sparse(1:Nx,1:Nx,[cc(2:end); cc(end-1)],Nx,Nx);
  
  Rxmu = (1/18/dx).*xMats.D3'*xMats.C3*B3*xMats.D3...
    + (1/144/dx).*xMats.D4'*xMats.C4*c*xMats.D4;
end
Mmux = xMats.D1int'*c*xMats.H*xMats.D1int + Rxmu;


% If boundaries are traction, then the formula reduces to this:
A = -xMats.Hinv*Mmux;


% Using the full formula, which is more easily generalized to different
% boundary conditions:
% muxBSx = c*xMats.BS;
% Dxxmu = xMats.Hinv*(-Mmux + muxBSx);
% 
% 
% % SAT penalty weight for traction BCs
alphaT = -1;
% 
% % add boundary condition terms
E0x = sparse(1,1,1,Nx,Nx);
ENx = sparse(Nx,Nx,1,Nx,Nx);
% 
% % traction
% AT = alphaT*xMats.Hinv*E0x*muxBSx;
% AB = alphaT*xMats.Hinv*ENx*muxBSx;

% A*u + (AT + AB)*gamma
AT = -alphaT*xMats.Hinv*E0x*c;
AB = alphaT*xMats.Hinv*ENx*c;
A2 = AT + AB;

% A = Dxxmu + AT + AB;
end

