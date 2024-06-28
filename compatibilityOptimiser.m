%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [outputSetIndexes,outputSetWeight] = compatibilityOptimiser(setIndexes,weight,setSize)
%compatibilityOptimiser finds the optimal combination by using maximising the inverse of the sum of the inverse of the weights of the sets

%input:
%   setIndexes          list of set vectors to be combined
%   weight              set weighting
%   setSize             number of elements to each set

%output:
%   output1             set combinations chosen by algorithm
%   output2             


if isscalar(weight)     % only one option
    outputSetIndexes = setIndexes{1};
    outputSetWeight = weight;
else                    %perform optimisation
    compatibleCombinationStorage = {};
    compatibleCombinationIndexes = {};
    compatibleCombinationWeight = {};

    tPC = numel(weight);%totalParameterCombinations 
    for a=1:tPC
        compatibleCombinationStorage{a} = [];
        compatibleCombinationWeight{a} = [];
        if a ==1
            combinationCount = nchoosek(tPC,1) ;
            combinations = nchoosek(1:tPC,1) ;
            for c=1:tPC
                compatibleCombinationStorage{a} = [ compatibleCombinationStorage{a} ; setIndexes{c} ] ;
                compatibleCombinationWeight{a} = [ compatibleCombinationWeight{a} ; weight(c) ] ;
            end
        else
            combinationCount = remainingCount*(nchoosek( tPC - a + 1 , 1 )) ;
            for c=1:remainingCount-1
                initialCombination = compatibleCombinationStorage{a-1}(c,:) ;
                initialWeight = compatibleCombinationWeight{a-1}(c);
                for d=c+1:tPC
                    if numel( unique([ initialCombination , setIndexes{d} ]) ) == numel ( [initialCombination , setIndexes{d} ] )
                        compatibleCombinationStorage{a} = [compatibleCombinationStorage{a} ; initialCombination , setIndexes{d}] ;
                        compatibleCombinationWeight{a} = [compatibleCombinationWeight{a} ; initialWeight + weight(d) ] ;
                        %save combination

                    else
                        %combination has a clash
                    end

                end
            end

        end
        remainingCount = numel(compatibleCombinationWeight{a}) ;
        if remainingCount == 0
            break
        end
    end

    maxWeight = 0;
    maxIndex = 0;
    maxCount = 0;
    for a=1:numel(compatibleCombinationWeight)
        [maxWeightLocal , maxIndexLocal] = max(compatibleCombinationWeight{a});
        if maxWeightLocal > maxWeight
            maxWeight = maxWeightLocal;
            maxIndex = maxIndexLocal;
            maxCount = a;
        end
    end

    outputSetIndexes = reshape( compatibleCombinationStorage{maxCount}(maxIndex,:) , [ maxCount , setSize ] ) ;
    outputSetWeight = maxWeight ;

end

end

