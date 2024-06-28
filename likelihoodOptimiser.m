%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [output] = likelihoodOptimiser(m,vB,dGrid,grid,precision,results,roundLimit)
%likelihoodOptimiser attempts to get the best possible likelihood function for the sum of parameters by optimising the parameter combinations. Calls recursiveSumFinder, likelihoodsArray, 
%compatibilityOptimiser, ArrayMultiplicity and likelihoodCombiner.

%input:
%   m               number of parameters involved in each measurement
%   vB              total number of parameters
%   dGrid           support grid spacing
%   precision       grid approximation precision
%   results         input results array
%   roundLimit      for time efficiency, influences algorithm choice. If possible combination count is below limit all treated at once, otherwise result count filtering is used.

%output:
%   output          optimised likelihood function

%further details:
%  By splitting the parameter estimation into the smallest combinations available and using the sets of combinations that have the most data we improve the limited data likelihood function creation
%  significantly compared to a brute force method of combining all at once.

if nargin<7 || isempty(roundLimit)
    roundLimit = 100;
end

total = nchoosek(vB,m);
multiple = nchoosek(vB-1,m-1);
possibilities = [];

for a=1:multiple
    if ~mod(multiple,a) && ~mod(total,a)
        possibilities = [possibilities;a,total/a,multiple/a];
    end
end

bestPossibility = possibilities(end,:); %number of sets    size of sets    multiple of theta'

numberOfSets = bestPossibility(1);
sizeOfSets = bestPossibility(2);
multipleOfTheta = bestPossibility(3);

parameterIndexList = zeros(numberOfSets,sizeOfSets);

%sort by results 5
[~,idx]=sort(results(5,:),'descend');

possibleParameterValues = flip(unique(results(5,:)));
possibleParameterCounts = zeros(numel(possibleParameterValues),1) ;
for a=1:numel(possibleParameterValues)
    possibleParameterCounts(a) = sum(results(5,:) == possibleParameterValues(a) );
end

parameterCountCategories= zeros(1,length(possibleParameterCounts));

%Defining the list of all of the parameters possible and padding with
%zeros.
allParametersUnpadded = nchoosek(1:vB,m);
allParameters = zeros(nchoosek(vB,m),vB);
for a=1:nchoosek(vB,m)
    for b=1:m
        allParameters( a , allParametersUnpadded(a,b) )= allParametersUnpadded(a,b);
    end
end

%adding a list to count the permutations in order. This will be used as a
%reference throughout. Will not clash with functions like sextor due to
%each having a unique identification number
temp = permute(1:nchoosek(vB,m),[2 1]);
allParameters = [allParameters , temp];
clear temp

%if expected possible combinations below upper limit perform one step algorithm
if nchoosek( nchoosek(vB,m) , sizeOfSets ) <= roundLimit
    parameters = cell(1);
    parameters{1} = allParameters;
    optimisationRounds = 1;
%else expected possible combinations too large so for time efficiency sort and perform calculations by result number
else
    parameters = cell(numel(possibleParameterCounts),1);

    %Sorting the results into sets of paramaters that correspond to the
    %ordered list of total result numbers
    for a=1:length(possibleParameterCounts)
        for b=1:nchoosek(vB,m)
            if results(5,b) == possibleParameterValues(a)
                parameters{a} = [parameters{a};allParameters(b,:)];
            end
        end

    end
    optimisationRounds = length(possibleParameterCounts) ;
end
totalCount = 0;
%maximisation step:

%declare the parameters that remain available
availableParameters = [];

choiceOuter = [] ;
choiceWeightOuter = [] ;

%loop over the different result counts. Begin with highest counts then
%add next set in each loop
for a=1:optimisationRounds
    %if final parameter has no results cut the end of the list of
    %combnations and break for loop
    if possibleParameterValues(a) == 0
        parameterIndexList(totalCount+1:end,:) = [];
        break
    end
    %add new set of parameters to the available parameters
    availableParameters = [availableParameters ; parameters{a}] ;
    total = multipleOfTheta*(1:vB);

    [availableParametersSize,~] = size(availableParameters);
    %declare arrays
    resultParametersStorage = cell( possibleParameterCounts(a) );
    indexToBeKept = zeros(possibleParameterCounts(a),sizeOfSets);
    %loop over the results in the newly added set. Any new result parameters
    %must include at least one parameter parameter from the new set.
    %Therefore we only need to search for new parameters that include at
    %least one of these newly considered parameters.

    options = cell(numel(possibleParameterCounts(a)),1);
    optionSize = zeros(numel(possibleParameterCounts(a)),1);
    optionWeight = zeros(numel(possibleParameterCounts(a)),1);
    optionInverseWeight = zeros(numel(possibleParameterCounts(a)),1);

    totalParameter = multipleOfTheta*(1:vB);
    repeat = 1;
    while repeat == 1
        [~,availableParametersTemp,repeat] = recursiveSumFinder( availableParameters(:,1:end-1) , totalParameter , sizeOfSets , availableParameters(:,end) ,1009);%max(24-3*vB,1) ); %1000);% max(24-3*vB,1) ) ;

        %Weight is the effective number of rounds if it was single parameter
        %estimation of theta (or ntheta ? where n=multipleOfTheta)
        possibleParameterCombinations = cell(length(availableParametersTemp),1);%zeros(numel(availableParametersTemp),sizeOfSets) ;
        weight = zeros(numel(availableParametersTemp),1);
        for b=1:numel(availableParametersTemp)
            for c=1:sizeOfSets
                possibleParameterCombinations{b}(c) = availableParametersTemp{b}(c,end) ;
                weight(b) = weight(b) + 1/results( 5, availableParametersTemp{b}(c,end) ) ;
            end
            weight(b) = 1/weight(b);
        end

        if ~isempty(weight)
            [ choice , choiceWeight ] = compatibilityOptimiser(possibleParameterCombinations,weight , sizeOfSets ) ;
            if isempty(choice)
                disp('choice empty')
            end
            for c=1:numel(choice)
                rowToRemove = find( availableParameters(:,end) == choice(c) ) ;
                availableParameters = setxor( availableParameters , availableParameters(rowToRemove, :) , 'rows' );
            end
            choiceOuter = [choiceOuter ; choice];
            choiceWeightOuter = [ choiceWeightOuter ; choiceWeight ];
        end

    end
end

%% Likelihood steps

llhOuter = zeros(numel(choiceWeightOuter) , precision ) ;
for a=1:numel(choiceWeightOuter)
    for b=1:sizeOfSets
        resultsTemp(:,b) = results( : , choiceOuter(a,b) ) ;
    end
    llhTemp = likelihoodsArray(sizeOfSets , grid , precision , resultsTemp) ;

    llhOuter(a,:) = arrayMultiplicity ( likelihoodCombiner(  llhTemp , sizeOfSets , precision )       , multipleOfTheta , precision , dGrid );
end

output = zeros(1,precision);
for a=1:precision
    output(a) = prod(llhOuter(:,a));
end

output = output/(sum(sort(output))*dGrid);

end



