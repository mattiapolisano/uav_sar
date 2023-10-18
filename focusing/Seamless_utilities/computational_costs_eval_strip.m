function [algo, varargout] = computational_costs_eval_strip(Sx,Sz, X, win, rho_x, lambda, OSF, t_ax, xmax,ymax)
%   MGP
%   UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


TDBP  = @(tau, N, M)      tau*M*N;

COSTO_FISSO =   @(tau,win,N,M)  tau/win*N*M*win;

COSTO_XY    =   @(Nim, N, M)   Nim*N*M;


Ntau=numel(Sx);
height= mean(Sz);
[Nx,Ny]= size(X);
Nr=numel(t_ax);
A_strip= lambda/2/rho_x * sqrt(ymax^2+xmax^2+height^2);

smin=sin(lambda/2/rho_x);
Nc=numel((-2/lambda*smin:1/(win*(lambda/4))/2/OSF:2/lambda*smin));
kmax=round(log2(Ntau/win));

for ii=1:kmax
    sl=min([1, sin(lambda/2/rho_x)+(1-sin(lambda/2/rho_x)*ii/kmax)]);
    Nloop=numel((-2/lambda*sl:1/min([A_strip, (win*lambda/4*2^ii)])/2/OSF:2/lambda*sl));
    costo_singolo_ciclo(ii)=Nloop*M*2^(kmax-ii);
end
costo_fisso=COSTO_FISSO(Ntau,win,Nc,M);

for ii=0:kmax
    costo_xy(ii+1)=COSTO_XY(2^(kmax-ii), Nx, Ny);
end


tdbpxy=TDBP(Ntau,Nx,Ny);


ffbprs=costo_fisso + sum(costo_singolo_ciclo);
ffbp1s=costo_fisso + costo_xy(1);



for ii=1:kmax
    temp(ii)=costo_fisso+sum(costo_singolo_ciclo(1:ii))+costo_xy(ii+1);
end

saf=[ffbp1s, temp];


figure; plot(0:kmax, saf, 'kd');
hold on; plot(0:kmax, (0:kmax)*0+ffbprs, 'b*')
hold on; plot(0:kmax, (0:kmax)*0+ffbp1s, 'r-o'), grid on
xlabel('hierarchical steps'); ylabel('Computational costs')
title('computational costs comparison')
legend('Seamless and Flexible', 'Full (r,s) FFBP', 'hybrid FFBP: single step')
drawnow


[safmin, safstop]=min(saf);

algo=[ffbprs tdbpxy,safmin ];

[cc, num_algo]=min(algo);

algo=num_algo;

if(num_algo==3)
    varargout{1}=safstop;
else 
    varargout{1}=0;
end





end