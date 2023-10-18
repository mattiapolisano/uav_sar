function [R,S] = Cart2Rs(X,Y,Z,centerx,centery,centerz)
%   MGP
%   [R,S] = Cart2Rs(X,Y,Z,centerx,centery,centerz)
%   Proietta un sistema di rifermiento
%   da un sistema cartesiano a uno RS
%   [R,S]   :  sistema di coordinate (R,S)
%   [X,Y,Z] :  sistema di coordinate cartesiano
%   [centerx, centery, centerz]  : centroidi delle coordinate del 
%                                  nuovo sistema di riferimento

R    =   sqrt((X-centerx).^2+(Y-centery).^2+(Z-centerz).^2);
S    =   (X-centerx)./R;
end