function [center_new, center_new_temp]= calcolo_nuovi_centri(center_old,ind_new)
% M.G.P.
% [center_new, center_new_temp] = calcolo_nuovi_centri(center_old,ind_new)
% Questa funzione ritorna i nuovi centri in due vettori
% center_new      :    1xN_new
% center_new_temp :    1xN_old
% center_old      :    vecchi centri
% ind_new         :    indici di centri da mediare


N_new       =   size(ind_new,2);    % numero di nuovi centri

center_new  =   zeros(1,N_new);         % inizializzo i nuovi centri

for ii=1:N_new      % ciclo lungo i nuovi centri
    
    % medio i centri a coppie
    center_new(ii)    =   mean(center_old(ind_new(:,ii)));  

end


% creo un vettore di dimensioni 1xN_old
% contenente i nuvi centri ripetuti


% inizializzo il vettore
center_new_temp    =    [];  

% genero il vettore di centri dim (2xN_new)
center_new_temp    =    [center_new; center_new];

% genero un vettore colonna N_oldx1
center_new_temp    =    center_new_temp(:);

end