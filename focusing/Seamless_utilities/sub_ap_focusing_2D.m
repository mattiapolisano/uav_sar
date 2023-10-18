function [stack, centerx, centery, centerz, sx_ax] = sub_ap_focusing_2D(D,t_ax,f0,Sx,Sy,Sz,win,psi_max,zfoc,OSF,varargin)
%   MGP
%   Questa funzione ritorna uno stack di immagini (range, sin) in banda base   
%   [stack, centerx, centery, centerz, sx_ax] = sub_ap_focusing(D,t_ax,f0,Sx,Sy,Sz,win,psi_max,zfoc,OSF)
%   senza pesi in x
%   [stack, centerx, centery, centerz, sx_ax] = sub_ap_focusing(D,t_ax,f0,Sx,Sy,Sz,win,psi_max,zfoc,OSF,kx0,DeltaK)
%   con i pesi in x

flag     =      ~isempty(varargin); % 1 se non è vuoto
                                    % 0 se è 

c        =      physconst('Lightspeed');     %velocità luce 

lambda   =      c/f0;               % lunghezza d'onda

Ntau     =      numel(Sx);          % # slow time samples

ind      =      finestratura(Ntau,win);      %indici sottoaperture

N        =      size(ind,2);        % numero di sottoaperture

r_ax     =      c.*t_ax./2;         % asse dei range

df       =      1/(win*lambda/4);   % campionamento in radianti

sin_max  =      sind(psi_max);      % max(s)

% asse s a bassa risoluzione
sx_ax_coarse   =    ...
    (-2/lambda*sin_max:df/OSF/2:2/lambda*sin_max)*lambda/2;

[S_coarse,R_coarse]   =   meshgrid(sx_ax_coarse,r_ax);

%inizializzo i centroidi delle sottoaperture

centerx_old   =    zeros(N,1);   %centroidi lungo x
centery_old   =    centerx_old;  % centroidi lungo y
centerz_old   =    centerx_old;  % centroidi lungo z

%inizializzo lo stack di immagini a bassa risoluzione

stack_old=zeros(size(R_coarse,1),size(R_coarse,2),N);


% (flag=1)
if (flag)     % se varargin non è vuoto 
              % (voglio usare i pesi in kx)

    kx0    =    varargin{1,1};       % centroide dei numeri d'onda
    DeltaK =    varargin{1,2};    % banda dei numeri d'onda

for n=1:N  % ciclo lungo le sottoaperture

    %calcolo i centri delle sottoaperture
    centerx_old(n)   =  mean(Sx(ind(:,n)));    %centro wrt x
    centery_old(n)   =  mean(Sy(ind(:,n)));    %centro wrt y
    centerz_old(n)   =  mean(Sz(ind(:,n)));    %centro wrt z
    
    % proietto il sistema (r,s) globale a bassa risoluzione 
    % in un cartesiano (X,Y,Z) a bassa risoluzione locale

    [X,Y,Z]    =   Rs2Cart(centerx_old(n),centery_old(n),centerz_old(n),zfoc,R_coarse,S_coarse);
    
    % Focalizzazione della sottoapertura
    temp       =   TDBP(D(:,ind(:,n)),X,Y,Z,f0,t_ax,Sx(ind(:,n)),Sy(ind(:,n)),Sz(ind(:,n)),kx0,DeltaK);

    % proiezione del cartesiano locale 
    % in un sistema (r,s) locale
    [RR,SS]    =   Cart2Rs(X,Y,Z,centerx_old(n),centery_old(n),centerz_old(n));

    %banda base e aggiornamento dello stack
    stack_old(:,:,n)   =    temp.*exp(-1j*4*pi/lambda.*RR);
    
end


else % flag = 0 ( non sto usando i pesi in kx)

for n=1:N     % ciclo lungo le sottoaperure

    %calcolo i centri delle sottoaperture
    centerx_old(n)   =  mean(Sx(ind(:,n)));    %centro wrt x
    centery_old(n)   =  mean(Sy(ind(:,n)));    %centro wrt y
    centerz_old(n)   =  mean(Sz(ind(:,n)));    %centro wrt z
    
    
    % proietto il sistema (r,s) globale a bassa risoluzione 
    % in un cartesiano (X,Y,Z) a bassa risoluzione locale
    [X,Y,Z]    =   Rs2Cart(centerx_old(n),centery_old(n),centerz_old(n),zfoc,R_coarse,S_coarse);
    
    % focalizzazione delle sottoaperture
    temp       =   TDBP(D(:,ind(:,n)),X,Y,Z,f0,t_ax,Sx(ind(:,n)),Sy(ind(:,n)),Sz(ind(:,n)));
    
    % proiezione del cartesiano locale in un
    % sistema (r,s) locale
    [RR,SS]    =    Cart2Rs(X,Y,Z,centerx_old(n),centery_old(n),centerz_old(n));

    % banda base e aggiornamento dello stack
    stack_old(:,:,n)   =    temp.*exp(-1j*4*pi/lambda.*RR);
    
end


end
% adesso ritorno le variabili di interesse

stack     =   stack_old;     % stack di immagini in banda base


% centri delle sottoaperture
centerx   =   centerx_old;   % rispetto a x
centery   =   centery_old;   % rispetto a y
centerz   =   centerz_old;   % rispetto a z


sx_ax     =   sx_ax_coarse;  % asse s a bassa risoluzione

end