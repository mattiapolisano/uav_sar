function I= TDBP_bist_unit(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz,varargin)
%   MGP
%   I= TDBP_bist_unit(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz)
%   caso bistatico, no pesi
%   I= TDBP_bist_unit(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz,k0x,DeltaK)
%   caso bistatico, pesi in kx


c=physconst('Lightspeed');
flag=~isempty(varargin);

if(flag)
    
    kx0=varargin{1};
    DeltaK=varargin{2};

    lambda=c/f0;
    
    Drx=X-Srx;
    Dry=Y-Sry;
    Drz=Z-Srz;
    R3r=sqrt(Drx.^2+Dry.^2+Drz.^2);
    R2r=sqrt(Drx.^2+Dy.^2);
    
    
    Dtx=X-Stx;
    Dty=Y-Sty;
    Dtz=Z-Stz;
    R3t=sqrt(Dtx.^2+Dty.^2+Dtz.^2);
    R3=R3r+R3t;
    
    Kx=1/lambda.*(R2r./R3r).*(Dx./R2r);
    W=rectpuls((Kx-kx0)/DeltaK);
    Tau=2.*R3./c;
    I=W.*interp1(t_ax,D,Tau).*exp(+1j*2*pi*f0*Tau);
else
    Drx=X-Srx;
    Dry=Y-Sry;
    Drz=Z-Srz;
    Rr=sqrt(Drx.^2+Dry.^2+Drz.^2);
    
    Dtx=X-Stx;
    Dty=Y-Sty;
    Dtz=Z-Stz;
    Rt=sqrt(Dtx.^2+Dty.^2+Dtz.^2);
    
    R=Rr+Rt;
    Tau=2.*R./c;
    I=interp1(t_ax,D,Tau).*exp(+1j*2*pi*f0*Tau);
end



end

