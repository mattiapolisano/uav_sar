function  SaveHDF5(Data,FolderPath)
%   MGP
%   SaveHDF5(Data,FolderPath)
%   Questa funzione salva variabili complesse in formato HDF5

filename  = getname(Data)+  ".h5";
Location  = FolderPath + "\"+ filename; 
% salvo la parte reale 

database  = '/real';

h5create(Location,database,size(Data), 'Datatype',class(Data))
h5write(Location,ds,real(Data))

database  = '/imag';

h5create(Location,database,size(Data), 'Datatype',class(Data))
h5write(Location,ds,imag(Data))


h5disp(Location)



end