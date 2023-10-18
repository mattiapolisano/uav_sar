function [X,Y,Z] = Rs2Cart(centerx,centery,centerz,zfoc,R,S)
%   MGP
%   [X,Y,Z] = Rs2Cart(centerx,centery,centerz,zfoc,R,S)
%   Proietta un sistema di riferimento
%   da coordinate RS a cartesiane 
%   [X,Y,Z] : coordinate cartesiane
%   [centerx, centery, centerz] : centri del sistema di riferimento 
%   [R,S] : sistema di riferimento (R,S)
%   zfoc  : quota del piano Z

X           =    R.*S+centerx;
Z           =    0*X+zfoc;
S(abs(S)>1) =    sign(S((abs(S)>1)))*1;
C           =    sqrt(1-S.^2);
Cos_theta   =    (Z-centerz)./(R.*C);

Cos_theta(abs(Cos_theta)>1)    =    1;

Sin_theta   =     sqrt(1-Cos_theta.^2);

Y           =     R.*C.*Sin_theta+centery;

end