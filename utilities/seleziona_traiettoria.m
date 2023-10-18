function [ind_start, ind_stop] = seleziona_traiettoria(f0,DeltaK,kx0,X,Y,Sx)
%   MGP
%   [ind_start, ind_stop] = seleziona_traiettoria(f0,DeltaK,kx0,X,Y,Sx)
%   Detailed explanation goes here
lambda     =  physconst('lightspeed')/f0;
rho_x      =  1/DeltaK;

psi        =  lambda/rho_x/2;

theta_hat  =  asin(kx0*lambda/2);

xmin       =  min(X(:));
xmax       =  max(X(:));
ymin       =  min(Y(:));
ymax       =  max(Y(:));


if (theta_hat==0)
    startx   = xmin  -  ymax*tan(psi/2);
    stopx    = xmax  -  ymax*tan(-psi/2);

elseif(theta_hat>=0)
    startx   = xmin  -  ymax*tan(theta_hat + psi/2);
    stopx    = xmax  -  ymin*tan(theta_hat - psi/2);

elseif(theta_hat<=0)

    startx   = xmin  -  ymin*tan(theta_hat + psi/2);
    stopx    = xmax  -  ymax*tan(theta_hat - psi/2);

end

[ ~ ,ind_start]  = min(abs(Sx-startx));
[ ~ ,ind_stop]   = min(abs(Sx-stopx));
end