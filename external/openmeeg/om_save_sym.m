function om_save_sym(data,filename,format)

% OM_LOAD_SYM   Load symmetric Matrix
%
%   Load symmetric Matrix
%
%   SYNTAX
%       [DATA] = OM_LOAD_SYM(FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' (default)
%
%   Created by Alexandre Gramfort on 2007-11-27.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.

if nargin < 3
    format = 'binary';
end

switch format
case 'binary'
    file = fopen(filename,'w');
    dims = size(data);
    fwrite(file,dims(1),'uint32','ieee-le');
    data = data(triu(ones(dims)>0));
    fwrite(file,data,'double','ieee-le');
    fclose(file);
case 'ascii'
    data = double(data);
    data = data(triu(ones(dim,dim)>0));
    save(filename,'data','-ASCII','-double','-v6')
otherwise
    error([me,' : Unknown file format'])
end
