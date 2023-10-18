function output = LoadHDF5(FileName)
%   MGP
%   output = LoadHDF5(FileName)  
%   Questa funzione legge file in formato .h5

database1  = '/real';
database2  = '/imag';

output=h5read(FileName,database1)+1j*h5read(FileName,database2);

end