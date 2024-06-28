%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [output] = arrayMultiplicity(inputArray,multiplicity,arrayPrecision,dGrid)
%arrayMultiplicity squeezes an array to repeat it over the same support multiplicity times, to be effective arrayPrecision >> multiplicity

%input:
%   inputArray          array to be repeated
%   multiplicity        how many repetitions in 
%   arrayPrecision      number of points initial array is defined over
%   dGrid               optional, allows for renormalisation of histograms

%output:
%   output              repeated array


if numel(inputArray) ~= arrayPrecision
    error('arrayMultiplicity input array precision')
end

temp = [];
for a=1:multiplicity
    temp = [temp,inputArray];
end
output = 0*inputArray;
for a=1:arrayPrecision
    output(a) = temp(multiplicity*a);
end

if ~isempty(dGrid)
    output = output/(sum(sort(output))*dGrid);
end

end

