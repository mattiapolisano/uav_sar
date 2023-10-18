function I = TDBP_unit(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz,varargin)
%   MGP
%   Questa funzione esegue la Time Domain Back Projection 
%   SOLO in caso monostatico
%   I = TDBP_unit(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz)
%   senza pesi in kx
%   I = TDBP_unit(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz,kx0,DeltaK)
%   con pesi in kx

flag    =    ~isempty(varargin); % 1 se non è vuoto
                                 % 0 se è vuoto

c       =    physconst('Lightspeed');   % velocità luce

if(flag)    % se sto usando i pesi 
    
    kx0    =     varargin{1};    % centrode kx
    DeltaK =     varargin{2};    % banda numeri d'onda
    
    lambda =     c/f0;           % lunghezza d'onda

    %Actual backprojection
    Dx    =     X-Sx;           % Delta X
    Dy    =     Y-Sy;           % Delta Y
    Dz    =     Z-Sz;           % Delta Z

    R3    =     sqrt(Dx.^2+Dy.^2+Dz.^2);   % Distanza in 3D
    R2    =     sqrt(Dx.^2+Dy.^2);         % Distanza in 2D

    Kx    =     2/lambda.*(R2./R3).*(Dx./R2); % frequenze d'onda
    
    %W    =     rectpuls((Kx-kx0)/DeltaK);    % pesi rettangolari
    
    W     =     exp(-0.5.*((Kx-kx0)./(DeltaK/2)).^2);  % pesi gaussiani
    
    Tau   =     2.*R3./c;       % Delay in 3D

    % interpolazione e banda passante
    I     =     W.*interp1(t_ax,D,Tau,'linear',0).*exp(+1j*2*pi*f0*Tau);
    
else                % non sto usando i pesi in x

    Dx   =   X-Sx;                      % delta X
    Dy   =   Y-Sy;                      % delta Y
    Dz   =   Z-Sz;                      % delta Z
    R    =   sqrt(Dx.^2+Dy.^2+Dz.^2);   % distanza in 3D
    Tau  =   2.*R./c;                   % delay in 3D

    % interpolazione e Banda Passante
    I    =   interp1(t_ax,D,Tau,'linear',0).*exp(+1j*2*pi*f0*Tau);
end

end

