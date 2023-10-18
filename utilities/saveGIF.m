function saveGIF(fig,filename,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    frame=getframe(fig);
    Im=frame2im(frame);
    
    [A,map] = rgb2ind(Im,256);
    if n == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",.5);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",.5);
    end

end