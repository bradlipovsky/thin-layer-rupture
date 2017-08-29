function [P, h] = NormalStress(t,M)
P = M.sigma;
h =  M.TideHeight * sin(2*pi*t/M.TidePeriod);
if M.Tides==2
    % TideH has the interpretation of a pressure head, 10 cm is about 1 kPa
    P = P + h;
end