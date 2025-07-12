% x^2 +(y/b)^2 <= (z+floor)^2/tan^2alpha

% b = 1; %before that we had 1.2
% tanAlpha = 2; %before that we had 1.1
% fl = 4;
% H = 2.1;

WGMs = 1.4;
WGMr = 1.8;
gamma = 1;

WGMh = gamma/WGMs;
cosPHI = 0.5*(WGMh + sqrt(2-WGMh^2));
DT = gamma*WGMr*sqrt(1 - cosPHI^2);
DB = 0.5*(gamma*WGMr*cosPHI - DT);
% DB = (gamma*WGMr*cosPHI);

H = DB + DT;
fl = WGMr*WGMs - DB;

b = 1; 
tanAlpha = WGMs; 
