function [stress,height] = TideStress(t,M)

g = 9.8 * 1e3; %mm/s^2   
rhow = 1; 

switch M.Tides
    case 1 % Longitudinal tides
        height = M.TideHeight * sin(2*pi*t/M.TidePeriod); % .* sin(2*pi*t/Ttide/14);
        stress = (M.rho*g*M.H^2 - rhow * g* (M.rho/rhow * M.H +height).^2 )/M.H/2;
    
  % case 2 is handeled in NormalStress.m
  
    case 3 % Longitudinal tides + geographic phase
        % Give the tide a spatially linear time lag in the arrival of the
        % tide.  Parameters result in a 4 hour time lag at sites 120km
        % apart.
        
        phase = M.TidePhase * linspace(0,1,M.N)';
        
        height = M.TideHeight * sin(2*pi*t/M.TidePeriod + phase); 
        stress = (M.rho*g*M.H^2 - rhow * g* (M.rho/rhow * M.H +height).^2 )/M.H/2;

    case 4 % Non zero longitudinal stress, but fixed zero tide height.
        ZeroTideStress = (M.rho*g*M.H^2 - rhow * g* (M.rho/rhow * M.H + 0).^2 )/M.H/2;
        stress = ZeroTideStress;
        
    case 5 % Ignore the background ice-static component.  Ignore Height^2 term because its tiny.
        height = M.TideHeight * sin(2*pi*t/M.TidePeriod);
        stress =  - M.rho *g * height;
        
    otherwise
        stress = 0;
        height=0;
end