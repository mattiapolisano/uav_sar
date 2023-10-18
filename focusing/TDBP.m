function I = TDBP(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz,varargin)
%   MGP
%   Questa funzione esegue la Time Domain Back Projection
%   I = TDBP(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz)
%   monostatico, no pesi
%   I = TDBP(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz, kx0,DeltaK)
%   monostatico, pesato in kx
%   I = TDBP(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz)
%   bistatico, no pesi
%   I = TDBP(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz, kx0,DeltaK)
%   bistatico, pesato

addpath(genpath('./TDBP_utilities'))

flag       =     ~isempty(varargin);   % 1 se non è vuoto
                                       % 0 se è vuoto

c          =     physconst('Lightspeed'); % velocità luce
lambda     =     c/f0;                    % lunghezza d'onda
Ntau       =     numel(Sx);               % # slow time

% inizializzo l'immagine focalizzata
I          =     0*X;


if(flag) % se non è vuoto
    % ci sono 3 casi:
    % - monostatico pesato in kx
    % - bistatico 
    % - bistatico pesato in kx

    if (numel(varargin)==2) % monostatico 
                            % pesato in kx

        kx0    = varargin{1}; % centroide kx
        DeltaK = varargin{2}; % banda kx 
        
        for ii=1:Ntau  % ciclo lungo l'apertura

            dat  =  D(:,ii);  % fast time signal
            % retroproiezione
            f    =  TDBP_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii),kx0,DeltaK);
            
            I    =  I+f; %somma immagini complesse
        end

    elseif(numel(varargin)==3)    % caso bistatico
                                  % non pesato in x
        % traiettorie trasmettitore
        Stx   =   varargin{1};  %wrt x
        Sty   =   varargin{2};  %wrt y
        Stz   =   varargin{3};  %wrt z

        for ii=1:Ntau  % ciclo lungo l'apertura

            dat   =  D(:,ii);   %fast time signal
            % retroproiezione
            f     =  TDBP_bist_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii),Stx(ii),Sty(ii),Stz(ii));
            
            I     =  I+f;       % somma immagini complesse
        end


    elseif(numel(varargin)==5)  % bistatico
                                % pesato in kx
        
        % traiettorie trasmettitore
        Stx    =   varargin{1};  %wrt x
        Sty    =   varargin{2};  %wrt y
        Stz    =   varargin{3};  %wrt z

        kx0    =   varargin{4};  % centroide pesi in kx
        DeltaK =   varargin{5};  % banda in kx 
        
        for ii=1:Ntau   % ciclo lungo l'apertura
            
            dat  =  D(:,ii);   %fast time signal
            % retroproiezione
            f    =  TDBP_bist_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii),Stx(ii),Sty(ii),Stz(ii),kx0,DeltaK);
            
            I    =   I+f;      % somma immagini complesse
        end

    else   % nessuno dei casi precendenti
        error('input error: wrong variable input')
    end


else       % monostatico
           % non pesato in kx

    for ii=1:Ntau    % ciclo lungo l'apertura
        
        dat   =   D(:,ii);  % fast time signal
        
        % retroproiezione
        f     =   TDBP_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii));
        
        I     =   I+f; % somma immagini complesse
    end
end


end

