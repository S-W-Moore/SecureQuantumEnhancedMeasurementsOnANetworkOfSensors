%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [outputSet,outputIndex,repeat] = recursiveSumFinder(inputArray,total,count,inputIndexes,outputLimit)
%recursiveSumFinder searches for the combinations of count input terms such that their sum equals the column number times total

%input:
%   inputArray         set of rows that could be combined
%   total               row to sum to
%   count               number of elements in combinations
%   inputIndexes        optional, variable indexing
%   outputLimit         for efficiency

%output:
%   outputSet             cell array of the combined sets
%   outputIndex           cell array of the indexes of the possible combinations

[numberOfRow, lengthOfRow] = size(inputArray);
repeat = 0;

outputSet = {};
outputIndex = {};
if isempty(inputIndexes)
    indexOuter = 1:1:numberOfRow;
else
    if  numberOfRow ~= length(inputIndexes)
        error('recursiveSumFinder2 input array and index lengths incompatible')
    end
    indexOuter = inputIndexes;
end
for a = 1:numberOfRow
    index = indexOuter;
    nestfun(inputArray(a,:),inputArray(a+1:end,:),index(a),index(a+1:end))
    if repeat == 1
        return
    end
end
    function nestfun(rowInitial,rowRemaining,indxInitial,indxRemaining)
        if ~all(sum(rowInitial,1)==total)
            for b = 1:numel(rowRemaining)/lengthOfRow
                nestfun([rowInitial;rowRemaining(b,:)],rowRemaining(b+1:end,:),[indxInitial;indxRemaining(b,:)],indxRemaining(b+1:end,:))
                if repeat == 1
                    return
                end
            end
        elseif any(sum(rowInitial,1)>total)
        elseif all(sum(rowInitial,1)==total) && numel(rowInitial) == count*lengthOfRow
            outputSet{end+1} = rowInitial;
            outputIndex{end+1} = indxInitial;
        end
        if length(outputIndex)>=outputLimit
            repeat = 1;
            return
        end
    end
end
