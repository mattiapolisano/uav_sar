function I= TDBP_bist_unit(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz,varargin)
%   MGP
%   Questa funzione esegue la Time Domain Back Projection
%   SOLO in caso bistatico
%   I= TDBP_bist_unit(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz)
%   caso bistatico, no pesi
%   I= TDBP_bist_unit(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz,k0x,DeltaK)
%   caso bistatico, pesi in kx

flag    =       ~isempty(varargin);   % 1 se non è vuoto
                                      % 0 se è vuoto

c       =       physconst('Lightspeed');   %  velocità luce


if(flag)    % pesi in kx
    
    kx0     =   varargin{1};       % centroide kx
    DeltaK  =   varargin{2};       % banda kx

    lambda  =   c/f0;              % lunghezza d'onda

    % calcolo le distanze scena-ricevitore
    Drx     =   X-Srx;             % Delta X
    Dry     =   Y-Sry;             % Delta Y
    Drz     =   Z-Srz;             % Delta Z
    R3r     =   sqrt(Drx.^2+Dry.^2+Drz.^2);   % Distanza 3D
    R2r     =   sqrt(Drx.^2+Dry.^2);           % Distanza 2D
    
    % calcolo le distanze scena-trasmettitore
    Dtx     =   X-Stx;             % Delta X
    Dty     =   Y-Sty;             % Delta Y
    Dtz     =   Z-Stz;             % Delta Z
    R3t     =   sqrt(Dtx.^2+Dty.^2+Dtz.^2);   % Distanza 3D


    Tau     =    (R3r+R3t)./c;     % ritardo totale
    % frequenza d'onda
    Kx      =    1/lambda.*(R2r./R3r).*(Drx./R2r);  

    W       =    rectpuls((Kx-kx0)/DeltaK);   % pesi freq. d'onda
    % focalizzazione immagine e banda passante
    I       =    W.*interp1(t_ax,D,Tau).*exp(+1j*2*pi*f0*Tau);


else      % senza pesi in kx 

    % calcolo le distanze scena-ricevitore

    Drx   =   X-Srx;     %Delta X
    Dry   =   Y-Sry;     %Delta Y
    Drz   =   Z-Srz;     %Delta Z
    Rr    =   sqrt(Drx.^2+Dry.^2+Drz.^2);      %distanza 3D

    % calcolo le distanze scena-trasmettitore    
    Dtx   =   X-Stx;     %Delta X
    Dty   =   Y-Sty;     %Delta Y
    Dtz   =   Z-Stz;     %Delta Z
    Rt    =   sqrt(Dtx.^2+Dty.^2+Dtz.^2);      %distanza 3D
    
    
    Tau   =   (Rr+Rt)./c;       % Delay
    % iterpolazione e banda passante
    I     =    interp1(t_ax,D,Tau).*exp(+1j*2*pi*f0*Tau);
end

end

