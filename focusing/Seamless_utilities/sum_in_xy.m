function Ixy= sum_in_xy(stack, X,Y,Z,centerx_old, centery_old, centerz_old, Sx,Sy,Sz, sx_old, r_ax, x_ax,y_ax)
%   MGP
%    Ixy= sum_in_xy(stack, X,Y,Z,centerx_old, centery_old, centerz_old, Sx,Sy,Sz, sx_old, r_ax, x_ax,y_ax)
%

N_im=size(stack,3);

    Ixy            =   0*X;
    zfoc           = mean(Z(:));
    for n=1:N_im

        [rr,ss]=Cart2Rs(centerx_old(n),centery_old(n),centerz_old(n), mean(Sx), mean(Sy), zfoc);
        sx_ap_temp=sx_old+ss;
        r_ax_temp= r_ax+rr;
        [rrr,sss]=meshgrid(r_ax_temp,sx_ap_temp);


        [xx,yy,zz]=Rs2Cart(mean(Sx),mean(Sy),mean(Sz),zfoc,rrr,sss);
        %         trovare il minimo
        [~, indxmin]=min(abs(x_ax-min(xx(:)-10)));
        [~, indymin]=min(abs(y_ax-min(yy(:))));

        % trovare il massimo
        [~, indxmax]=min(abs(x_ax-max(xx(:)+10)));
        [~, indymax]=min(abs(y_ax-max(yy(:))));
        % selezionarli per farci l'interpolazione
        indx=(indxmin:indxmax);
        indy=(indymin:indymax);


        [Rxy,Sxy]  =   ...
            Cart2Rs(X(indy,indx),Y(indy,indx),Z(indy,indx),centerx_old(n),centery_old(n),centerz_old(n));
        temp       =    interp2(sx_old,r_ax,stack(:,:,n),Sxy,Rxy,'cubic',0).*exp(+1j*4*pi/lambda.*Rxy);

        % somma
        Ixy(indy,indx)        =    Ixy(indy,indx) +temp;

    end

end