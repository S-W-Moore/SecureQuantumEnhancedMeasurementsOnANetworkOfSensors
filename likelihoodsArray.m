%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285

function array = likelihoodsArray(llhNumber,grid,precision,results)
%likelihoodsArray ouputs an array as a set of grid approximation to a set of likelihood functions

%input:
%   llhNumber           how many likelihoods to be computed
%   grid                grid of likelihood funciton support
%   precision           grid points
%   results             measurement result count

array = zeros(llhNumber,precision);
dGrid = grid(2)-grid(1);
for a=1:llhNumber
    %multiplier allows the analysis of larger data.
    multiplier = 1;
    if results(5,a) > 500
        multiplier = ceil( results(5,a) / 500 );
    end
    tempResults = double(results(:,a)) / multiplier;
    array(a,:) = (.25*(1+cos(grid))).^tempResults(1) .* (.25*(1+cos(grid+pi/2))).^tempResults(2) .* (.25*(1+cos(grid+pi))).^tempResults(3) .* (.25*(1+cos(grid+3*pi/2))).^tempResults(4);
    array(a,:) = array(a,:)/sum(sort(array(a,:))*dGrid);
    if multiplier>1
        array(a,:) = array(a,:).^multiplier;
        array(a,:) = array(a,:)/sum(sort(array(a,:))*dGrid);
    end
end

end