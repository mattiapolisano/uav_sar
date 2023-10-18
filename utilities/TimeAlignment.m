function t_ax_new = TimeAlignment(t_ax,D)
%   MGP
%    t_ax_new = TimeAlignement(t_ax,D,t0)
%   Detailed explanation goes here


[t,kt]=max(D);

figure; plot(t_ax(kt)); xlabel('slow time bins'); ylabel('first arrival time');
title('first arrival time')

t_average=mean(t_ax(kt));

t_ax_new=t_ax-t_average;

end