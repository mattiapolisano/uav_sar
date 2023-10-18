function [ind,Sx_proc,Sy_proc,Sz_proc,D_proc] =NextPower2Selection(ind_start, ind_stop,lambda,Sx,Sy,Sz,D)

%   [ind,Sx_proc,Sy_proc,Sz_proc,D_proc] =NextPower2Selection(ind_start, ind_stop,lambda,Sx,Sy,Sz,D)
%   Detailed explanation goes here

indind       =  [ind_start, ind_stop];
ind          =  (min(indind): max(indind));
indmean      =  2^floor(log2(mean(ind)));
N_proc       =  ind_stop-ind_start;
N_pow_2      =  2^nextpow2(N_proc);


Sx_proc      =  Sx(indmean-N_pow_2/2+1:indmean+N_pow_2/2);

Sy_proc      =  Sy(indmean-N_pow_2/2+1:indmean+N_pow_2/2);

Sz_proc      =  Sz(indmean-N_pow_2/2+1:indmean+N_pow_2/2);

D_proc       = D(:,indmean-N_pow_2/2+1:indmean+N_pow_2/2);
end