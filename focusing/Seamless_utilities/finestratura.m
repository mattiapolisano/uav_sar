function indici = finestratura(N,l,varargin)
%   MGP
%   indici = finestratura(N,l)
%   crea un set di indici di N/l finestre
%   N  :  numero di elementi DA FINESTRARE
%   l  :  lunghezza della finestra
%   (Da finire )
%   indici = finestratura(N,l, overlap_samples)
%   crea un set di indici con overlap tra le finestre.
%   overlap_samples  : è espresso in numero di campioni


flag       =   ~isempty(varargin);  % 1 se non è vuoto
% 0 se è

if (flag)  % se voglio usare l'overlapp
    overlap_samples    =   l-varargin{1};
    % finestratura con overlap
    indici             =   [];    % inizializzo gli indici

    jj                 =   0;     % iteratore

    for ii=l:overlap_samples:N    % ciclo lungo le aperture

        jj             =   jj+1;  % update iteratore

        indici(:,jj)   =   (ii-l+1:ii);   % aggiorno gli indici

    end

else    % non sto usando overlapp


    L       =      N/l;                  % numero di finestre
    indici  =      zeros(l,L);           % inizializzo gli indici

    for n=1:L  % itero lungo il numero di finestre

        % aggiorno gli indici
        indici(:,n)     =     ((n-1)*l+1:(n-1)*l+l);

    end

end
end

