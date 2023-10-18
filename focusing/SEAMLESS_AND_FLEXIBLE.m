function [Ixy, asse] = SEAMLESS_AND_FLEXIBLE(D,t_ax,f0,Sx,Sy,Sz,win,rho_x,zfoc,OSF,X,Y,Z,type, varargin)
%   MGP
%   [Irs, s_ax] = SEAMLESS_AND_FLEXIBLE(D,t_ax,f0,Sx,Sy,Sz,win,rho_x,zfoc,OSF,X,Y,Z,type)
%   performs the focusing of a area defined in [X,Y,Z] and returns a
%   focused image in (r,s). The focusing is perfomed without spectral
%   weights
%   [Ixy, asse] = SEAMLESS_AND_FLEXIBLE(D,t_ax,f0,Sx,Sy,Sz,win,rho_x,zfoc,OSF,X,Y,Z,type, kx0, DeltaK)
%   returns the image in (x,y) ot (r,s)according to faster focusing algorithm.

addpath(genpath('./Seamless_utilities'))


if(numel(varargin)==2)
    kx0=varargin{1};
    DeltaK=varargin{2};
else
    error('wrong input arguements')
end
c=physconst('Lightspeed');
lambda=c/f0;
psi_max=90;

if(strcmp(type, 'slitta'))

    flagxy=false;
    [Ixy, asse] = seamless_and_flexible_slitta(D,t_ax,f0,Sx,Sy,Sz,win,psi_max,zfoc,OSF, X,Y,Z, flagxy);

elseif(strcmp(type, 'stripmap'))

    [algo, stepout] = computational_costs_eval_strip(Sx,Sz, X, win, rho_x, lambda, OSF, t_ax, xmax,ymax);
    [Ixy, asse] = seamless_and_flexible_stripmap(D,t_ax,f0,Sx,Sy,Sz,win,rho_x,zfoc,OSF,X,Y,Z,kx0,DeltaK,algo,stepout);

else
    error('Wrong type: or  "slitta" or "stripmap" ')
end

end