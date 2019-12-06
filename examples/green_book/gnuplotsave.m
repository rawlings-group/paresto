function gnuplotsave(datname, data, fmt)
% gnuplotsave(datname, data, [fmt='%g'])
%
% Saves the fields of struct 'data' to the given filename. This script uses a
% format that allows gnuplot to use named blocks.
%
% fmt is the format string to use for numeric outputs. The efault '%g' includes
% 6 significant figures. Use, e.g., '%.10g' for 10 significant figures.
narginchk(2, 3);
if nargin() < 3
    fmt = '%g';
end
if isstruct(data)
    method = 'struct';
elseif ismatrix(data)
    method = 'matrix';
else
    error('Unknown data type: %s', type(data));
end
if exist(datname, 'file')
    delete(datname);
end
switch method
case 'struct'
    savestruct(datname, data, '', fmt);
case 'matrix'
    save('-ascii', datname, 'data');
end

end%function

function savestruct(datname, data, prefix, numberformat)
% savestruct(datname, data, prefix, numberformat)
%
% Saves struct fileds to a data file. Recursive for nested structs.
narginchk(4, 4);
fields = fieldnames(data);
for i = 1:length(fields)
    f = fields{i};
    val = data.(f);
    valname = [prefix, f];
    if isstruct(val)
        savestruct(datname, val, [valname, '.'], numberformat);
    elseif isnumeric(val) && ismatrix(val)
        dat = fopen(datname, 'a');
        fprintf(dat, '# [%s]\n', valname);
        fmt = [repmat([numberformat, ' '], 1, size(val, 2)), '\n'];
        fprintf(dat, fmt, double(val)');
        fprintf(dat, '\n\n');
        fclose(dat);
    else
        error('Field %s.%s is invalid. Must be matrix or struct!', prefix, f);
    end
end

end%function

