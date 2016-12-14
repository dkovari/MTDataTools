function c = cell_sprintf(format,data)
% sprintf with cellstr array output
%   format: c-style format string
%   data: array of data to use in formatted string


%determine number of outputs in string
nOut = numel(regexp(format,'%[-+0-9cdeEfgGosuxX]'));

flmat = false;
if ismatrix(data) && size(data,2)==nOut
    data = data';
    flmat = true;
end

reshape(data,nOut,[]);
c = cell(1,size(data,2));
for n=1:size(data,2)
    c{n} = sprintf(format,data(:,n));
end

if flmat
    c = c';
end