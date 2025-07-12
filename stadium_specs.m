function [s,H,fl,l] = stadium_specs(tau)
    % % H = tau;
    H = 1;
    s = tan(tau);
    % s = tan(pi/3);
    topX = 2*sqrt(2);
    % topX = 2.2;
    fl = 0.5*topX*s - H;
    l = 0.6;

