function c = cell_sprintf(format,data)
% sprintf with cellstr array output
%   format: c-style format string
%   data: array of data to use in formatted string

c = {};
for dataElement = data
    c{end+1} = sprintf(format,dataElement);
end
end