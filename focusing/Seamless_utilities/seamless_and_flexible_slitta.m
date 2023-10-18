function [I, asse] = seamless_and_flexible_slitta(D,t_ax,f0,Sx,Sy,Sz,win,psi_max,zfoc,OSF, X,Y,Z, flagxy)
%   MGP
%   I = seamless_and_flexible_slitta()
%   Detailed explanation goes here

l=2;
c=physconst('lightspeed');
lambda=c/f0;
r_ax=t_ax*c/2;
x_ax=X(1,:);
y_ax=Y(:,1);


[stack_old, centerx_old, centery_old, centerz_old, sx_ax_coarse] = sub_ap_focusing_2D(D,t_ax,f0,Sx,Sy,Sz,win,psi_max,zfoc,OSF);


A_win   =    win*lambda/4;   %lunghezza della sottoapertura
N=size(stack_old,3);
Ntau=numel(Sx);
ind=finestratura(Ntau,win);
% inizializzo le variabili del loop
A_old   =    A_win;           % lunghezza ap. al passo n-1
N_old   =    N;               % numero sotto ap. al passo n-1
sx_old  =    sx_ax_coarse;    % asse s al passo n-1
ind_old =    ind;             % indici delle sottoaperture
% al passo n-1
stack   =    stack_old;       % stack di immagini in Banda Base
% al passo n-1

maxstep = log2(Ntau/win);
step    = 0;

while(size(ind_old,2)>1) % ciclo finché non rimane solo un'immagine

    %calcolo i nuovi assi e le nuove coppie
    step    = step+1;
    % la lunghezza delle sottoaperture raddoppia (*l)
    A_new   =    l*A_old;
    sloop   =  1;
    % calcolo gli indici da sommare a coppie
    ind_new =    finestratura(N_old,l);
    % numero delle nuove sottoaperture
    N_new   =    size(ind_new,2);

    df      =    1/ min([  A_new]);        % campionamento in radianti

    % nuovo asse s
    sx_new  =    (-2/lambda*sloop:df/2/OSF:2/lambda*sloop)*lambda/2;

    [S_new,R_new]    =    meshgrid(sx_new,r_ax);

    % calcolo i nuovi centri delle aperture
    % creo un vettore con i nuovi centri
    % dimensionati giusti
    [centerx_new,centerx_new_temp]=...
        calcolo_nuovi_centri(centerx_old, ind_new);

    [centery_new,centery_new_temp]=...
        calcolo_nuovi_centri(centery_old, ind_new);

    [centerz_new,centerz_new_temp]=...
        calcolo_nuovi_centri(centerz_old, ind_new);


    % inizializzo lo stack di immagini
    % interpolate sulla griglia più fitta
    temp_stack=zeros(size(R_new,1),size(R_new,2),N_old);

    % interpolo e banda passante
    for n=1:N_old

        %proietto il sistema (r,s) in un
        %sistema cartesiano comune a coppie
        [X_new,Y_new,Z_new]   =...
            Rs2Cart(centerx_new_temp(n),centery_new_temp(n),centerz_new_temp(n),zfoc,R_new,S_new);
        %proietto il sistema cartesiano comune a coppie
        %su un sistema (r,s) locale della singola sottoapertura
        [R,S]   = ...
            Cart2Rs(X_new,Y_new,Z_new,centerx_old(n),centery_old(n),centerz_old(n));


        %interpolazione su griglia più fitta e banda passante
        temp_stack(:,:,n)=interp2(sx_old,r_ax,stack(:,:,n),S,R,'cubic',0).*exp(+1j*4*pi/lambda.*R);


    end

    %inizializzo lo stack
    stack_new=zeros(size(R_new,1),size(R_new,2),N_new);

    % somma a coppie e Banda Base
    for n=1:N_new

        %proietto il sistema (r,s) in un
        %sistema cartesiano comune a coppie

        [X_new,Y_new,Z_new]=...
            Rs2Cart(centerx_new(n),centery_new(n),centerz_new(n),zfoc,R_new,S_new);

        %proietto il sistema cartesiano comune a coppie
        %su un sistema (r,s) comune a coppie
        [RR,S] = ...
            Cart2Rs(X_new,Y_new,Z_new,centerx_new(n),centery_new(n),centerz_new(n));

        % somma e Banda Base
        stack_new(:,:,n)=...
            sum(temp_stack(:,:,ind_new(:,n)),3).*exp(-1j*4*pi/lambda.*RR);

    end


    % loop variables
    stack   =   stack_new;   % stack di immagini
    A_old   =   A_new;       % lunghezza sottoapertura
    N_old   =   N_new;       % numero di sottoaperture

    centerx_old   =   centerx_new;   %centroidi delle sottoaperture
    centery_old   =   centery_new;   %centroidi delle sottoaperture
    centerz_old   =   centerz_new;   %centroidi delle sottoaperture

    sx_old   =   sx_new;   % asse s
    ind_old  =   ind_new;  % indici delle sottoaperture

end

if(flagxy)

    I= sum_in_xy(stack, X,Y,Z,centerx_old, centery_old, centerz_old, Sx,Sy,Sz, sx_old, r_ax, x_ax,y_ax);
    asse=x_ax;

else

    I=stack.*exp(+1j*4*pi/lambda.*R_new);
    asse=sx_old;

end

end