function I = TDBP(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz,varargin)
%   MGP
%   I = TDBP(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz)
%   monostatico, no pesi
%   I = TDBP(D,X,Y,Z,f0,t_ax,Sx,Sy,Sz, kx0,DeltaK)
%   monostatico, pesato in kx
%   I = TDBP(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz)
%   bistatico, no pesi
%   I = TDBP(D,X,Y,Z,f0,t_ax,Srx,Sry,Srz,Stx,Sty,Stz, kx0,DeltaK)
%   bistatico, pesato
c=physconst('Lightspeed');
lambda=c/f0;
Ntau=numel(Sx);

I=0*X;

flag=~isempty(varargin);

if(flag)
    if (numel(varargin)==2)
        kx0=varargin{1};
        DeltaK=varargin{2};
        %caso2
        for ii=1:Ntau
            dat=D(:,ii);
            f=TDBP_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii),kx0,DeltaK);
            I=I+f;
        end
    elseif(numel(varargin)==3)
        %caso 3
        Stx=varargin{1};
        Sty=varargin{2};
        Stz=varargin{3};
        for ii=1:Ntau
            dat=D(:,ii);
            f=TDBP_bist_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii),Stx(ii),Sty(ii),Stz(ii));
            I=I+f;
        end
    elseif(numel(varargin)==5)
        %caso 4
        Stx=varargin{1};
        Sty=varargin{2};
        Stz=varargin{3};
        kx0=varargin{4};
        DeltaK=varargin{5};
        
        for ii=1:Ntau
            dat=D(:,ii);
            f=TDBP_bist_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii),Stx(ii),Sty(ii),Stz(ii),kx0,DeltaK);
            I=I+f;
        end
        
        
    else
        error('input error')
    end
else
    %caso 1
    for ii=1:Ntau
        dat=D(:,ii);
        f=TDBP_unit(dat,X,Y,Z,f0,t_ax,Sx(ii),Sy(ii),Sz(ii));
        I=I+f;
    end
end


end

