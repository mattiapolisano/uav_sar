function indici = finestratura(N,l,varargin)
%   MGP
%   indici = finestratura(N,l)
%   crea un set di indici, N/l finestre
%   N numero di elementi DA FINESTRARE
%   l lunghezza della finestra
%   indici = finestratura(N,l, overlap_samples)
%   crea un set di indici con overlap tra le finestre.
%   overlap_samples è espresso in numero di campioni


flag=~isempty(varargin);
if (flag)
    overlap_samples=varargin{1};
    % finestratura con overlap
    indici=[];
    jj=0;
    
    for ii=l:overlap_samples:N
        
        jj=jj+1;
        indici(:,jj)=(ii-l+1:ii);
        
    end
    
else
    
    
    L=N/l;
    indici=zeros(l,L);
    
    for(n=1:L)
        indici(:,n)=((n-1)*l+1:(n-1)*l+l);
        
    end
    
end




end

