function [v,a,ter] = calc_ddm(pc,vrt,mrt)

% calculates diffusion parameters from accuracy and response time distrib.
% adpated from (http://www.ejwagenmakers.com/2007/EZ.pdf)
% http://www.ejwagenmakers.com/EZ.html

% inputs
% pc = percent correct (min=0,max=1)
% var = var(rt) (response time in ms)
% mrt = mean(rt)

% output
% v = drift rate
% a = boundary separation
% ter = nondecision time

% transform to
% seconds squared
vrt = (vrt/10^6);

% transform 
% to seconds 
mrt = (mrt/10^3);

if length(pc)>1
    pc = mean(pc);
end

s = .1; % arbitrary scaling parameter  
s2 = s^2; % violitility parameter 

% edge correction 
if pc == 0 || .5 || 1
    pc=pc-(rand/100); 
end

% calculate 
% diffusion 
% parameters 
L = qlogit(pc);
x = L*(L*pc^2-L*pc+pc-.5)/vrt; 
v = sign(pc-.5)*s*(x)^(1/4); % drift rate 
a = s2*L/v; % response boundary 
y = -v*a/s2;
mdt = (a/(2*v))*(1-exp(y))/(1+exp(y));
ter = mrt-mdt; % non-decision time 

end