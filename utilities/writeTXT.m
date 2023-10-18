function  writeTXT(experiment_folder,radar_parameters)
%   writeTXT(experiment_folder,radar_parameters)
%   Detailed explanation goes here

names=fieldnames(radar_parameters);
values=struct2cell(radar_parameters);
names(8:9)=[];
values(8:9)=[];
fid=fopen(experiment_folder+'\radar_parameters.txt','wt');
for ii=1:numel(names)
    fprintf(fid, '%s : %0.001f \n', cell2mat(names(ii)), cell2mat(values(ii)));
end
fclose(fid);
end