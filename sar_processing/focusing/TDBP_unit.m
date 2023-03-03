function I = TDBP_unit(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz,varargin)
%   MGP
%   I = TDBP_unit(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz)
%   caso monostatico, no pesi
%   I = TDBP_unit(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz,kx0,DeltaK)
%   caso monostatico, pesi in x

c=physconst('Lightspeed');
flag=~isempty(varargin);

if(flag)
    
    kx0=varargin{1};
    DeltaK=varargin{2};
    
    lambda=c/f0;
    Dx=X-Sx;
    Dy=Y-Sy;
    Dz=Z-Sz;
    R3=sqrt(Dx.^2+Dy.^2+Dz.^2);
    R2=sqrt(Dx.^2+Dy.^2);
    Kx=2/lambda.*(R2./R3).*(Dx./R2);
    W=rectpuls((Kx-kx0)/DeltaK);
    Tau=2.*R3./c;
    I=W.*interp1(t_ax,D,Tau,'linear',0).*exp(+1j*2*pi*f0*Tau);
else
    Dx=X-Sx;
    Dy=Y-Sy;
    Dz=Z-Sz;
    R=sqrt(Dx.^2+Dy.^2+Dz.^2);
    Tau=2.*R./c;
    I=interp1(t_ax,D,Tau,'linear',0).*exp(+1j*2*pi*f0*Tau);
end


end

