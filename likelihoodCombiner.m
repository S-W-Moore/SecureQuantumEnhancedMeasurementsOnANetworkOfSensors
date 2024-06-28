%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [output] = likelihoodCombiner(likelihood,likelihoodNumber,likelihoodPrecision)
%likelihoodCombiner uses fast fourier transforms to create a likelihood function for a sum of parmaeters

%input:
%   likelihood          array of likelihoods to be combined
%   likelihoodNumber    number of likelihoods to be combined
%   likelihoodPrecision array precision of likelihoods

%output:
%   output              single likelihood for the sum of the parameters making each row of 'likelihood' input

%further information:
%   The use of fast fourier transforms allows us to perform this calculation at approximately NlogN rather than N^2 by brute force where N=likelihoodPrecision

if size(likelihood) ~= [likelihoodNumber,likelihoodPrecision]
    error('likelihoodCombiner likelihood precision incompatible with other inputs')
end
transformedLikelihoods = fft( likelihood(1,:) ) ;
if likelihoodNumber>1
    for a=2:likelihoodNumber
        %in fourier space the convolution of the likelihoods is the product of their fourier transforms
        transformedLikelihoods = transformedLikelihoods .* fft(likelihood(a,:));
    end
end
%invert fourier transforms
output = ifft(transformedLikelihoods);
%remove some numerical noise due to using fft and ifft
output (output<0)= 0;

end

