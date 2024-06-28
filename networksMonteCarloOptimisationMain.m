%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285

function networksMonteCarloOptimisationMain(nMax,EveLimit,attackTypeTotal,vBMin,vBMax,vBStep,precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue,outputUserPath)
%networksMonteCarloOptimisationMain runs the F and S optimisation protocol for separable initial states, entangled initial states and hydbrid.

%input:
%   nMax                        maximum rounds
%   EveLimit                    security limit placed on Eve
%   precisionRoot               data points in likelihood function grid approximations is 2^precisionRoot
%   attackTypeTotal             scalar or array, 0 Alice only wihtout attack, 1-4 different attacks by Eve on every round
%           1  -   replace with separable state
%           2  -   measure all in random basis send separable state
%           3  -   replace with entangled state
%           4  -   measure as if entangled then replace with entangled state
%   vBmin, vBMax, vBSte         minimum, maximum and steps in vB values
%   steps1                      initial grid size for both dimensions
%   steps2                      number of repetitions of new grid step
%   repetitionTrueValue         number of different true values tested
%   repetitionEachValue         number of repetitions for each true value
%   ouptutUserpath              destination for output files

%ouput:
%   saves files

if vBMax>127
    error('networksMonteCarloOptimisationMain optimisation requires vB to be less than int8')
end

%define output userpath
userpath(outputUserPath)

%simulation parameters
Smin1 = 0;
Smax1 = 1;
Fmin1 = 0;
Fmax1 = 1;
vSArray1 = linspace(Smin1,Smax1,steps1);
vFArray1 = linspace(Fmin1,Fmax1,steps1);
precision = 2^precisionRoot;
dPhi = 2*pi/precision;
Phi = linspace(dPhi, 2*pi, precision);

%reset random number stream
rng('shuffle','twister')

%creation of parameter values
vBArray = vBMin:vBStep:vBMax;
vBTotal = length(vBArray);
trueValueCell = cell(vBTotal,1);
for vBIndex=1:vBTotal
    trueValueCell{vBIndex} = 2*pi*rand( repetitionTrueValue , vBArray(vBIndex) );
end


%declaration of Round 1 storage variables
attackTypeCount = length(attackTypeTotal);
AliceFI = zeros(vBTotal,steps1,steps1);
AliceMSE = zeros(vBTotal,steps1,steps1);
EveMSE = cell(attackTypeCount,1);
roundCountAverage = cell(length(attackTypeTotal),1);
for attackTypeIndex = 1:length(attackTypeTotal)
    EveMSE{attackTypeIndex} = zeros(vBTotal,steps1,steps1);
    roundCountAverage{attackTypeIndex} = zeros(vBTotal,steps1,steps1);
    undetectionsAverage{attackTypeIndex} = zeros(vBTotal,steps1,steps1);
end

%declaration of Round 2 storage variables
AliceFI2 = cell(attackTypeCount,1);
AliceMSE2 = cell(attackTypeCount,1);
EveMSE2 = cell(attackTypeCount,1);
roundCountAverage2 = cell(attackTypeCount,1);
estimationPoints2 = cell(attackTypeCount,1);
for attackTypeIndex = 1:length(attackTypeTotal)
    AliceFI2{attackTypeIndex} = cell(vBTotal,1);
    AliceMSE2{attackTypeIndex} = cell(vBTotal,1);
    EveMSE2{attackTypeIndex} = cell(vBTotal,1);
    roundCountAverage2{attackTypeIndex} = cell(vBTotal,1);
    undetectionsAverage2{attackTypeIndex} = cell(vBTotal,1);
    estimationPoints2{attackTypeIndex} = cell(vBTotal,1);
end

outerMaxM = cell(attackTypeCount,1);
outerMaxE = cell(attackTypeCount,1);
outerMaxS = cell(attackTypeCount,1);
for attackTypeIndex = 1:attackTypeCount
    outerMaxM{attackTypeIndex} = zeros(vBTotal,7);
    outerMaxE{attackTypeIndex} = zeros(vBTotal,7);
    outerMaxS{attackTypeIndex} = zeros(vBTotal,7);
end

%begin time for tracking progress
tic

for vBIndex=1:vBTotal
    %%  First round, initial grid
    vB = vBArray(vBIndex);
    toc
    tic
    disp(['First round, vB = ',num2str(vB)])
    if vB == 1              %Avoid repetitions for vB = 1 which don't make sense
        steps1Temp = 1;
    else
        steps1Temp = steps1;
    end
    attackTypeLocal = [0,attackTypeTotal];
    for vSIndex=1:steps1Temp
        vS = vSArray1(vSIndex);
        for vFIndex=1:steps1
            vF = vFArray1(vFIndex);
            [AliceMSE(vBIndex,vSIndex,vFIndex),EveMSETemp,roundCountAverageTemp,AliceFI(vBIndex,vSIndex,vFIndex),undetectionsAverageTemp] = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,nMax,trueValueCell{vBIndex},precision,vB,vS,vF,attackTypeLocal) ;
            for attackTypeIndex = 1:length(attackTypeTotal)
                EveMSE{attackTypeIndex}(vBIndex,vSIndex,vFIndex) = EveMSETemp(attackTypeIndex);
                roundCountAverage{attackTypeIndex}(vBIndex,vSIndex,vFIndex) = roundCountAverageTemp(attackTypeIndex);
                undetectionsAverage{attackTypeIndex}(vBIndex,vSIndex,vFIndex) = undetectionsAverageTemp(attackTypeIndex);
            end
        end
    end

    %%  Second round, optimisation algorithm

    for attackTypeLocalIndex = 1:length(attackTypeTotal)
        attackTypeLocal = [0,attackTypeTotal(attackTypeLocalIndex)];
        toc
        tic
        disp(['Second round, vB = ',num2str(vB) , '; attackType = ',num2str(attackTypeLocal(2))])

        vFLength = steps1;
        if vB == 1
            vSLength = 1;
        else
            vSLength = steps1;
        end

        for vSIndex=1:steps1
            for vFIndex= 1:steps1
                AliceFI2{attackTypeLocalIndex}{vBIndex} = [ AliceFI2{attackTypeLocalIndex}{vBIndex} ; AliceFI(vBIndex,vSIndex,vFIndex) ] ;
                AliceMSE2{attackTypeLocalIndex}{vBIndex} = [ AliceMSE2{attackTypeLocalIndex}{vBIndex} ; AliceMSE(vBIndex,vSIndex,vFIndex) ] ;
                EveMSE2{attackTypeLocalIndex}{vBIndex} = [ EveMSE2{attackTypeLocalIndex}{vBIndex} ; EveMSE{attackTypeLocalIndex}(vBIndex,vSIndex,vFIndex) ] ;
                roundCountAverage2{attackTypeLocalIndex}{vBIndex} = [ roundCountAverage2{attackTypeLocalIndex}{vBIndex} ; roundCountAverage{attackTypeLocalIndex}(vBIndex,vSIndex,vFIndex) ] ;
                undetectionsAverage2{attackTypeLocalIndex}{vBIndex} = [ undetectionsAverage2{attackTypeLocalIndex}{vBIndex} ; undetectionsAverage{attackTypeLocalIndex}(vBIndex,vSIndex,vFIndex) ] ;
                estimationPoints2{attackTypeLocalIndex}{vBIndex} = [ estimationPoints2{attackTypeLocalIndex}{vBIndex} ; vSArray1(vSIndex) , vFArray1(vFIndex) ] ; 
            end
        end

        totalCount = steps1^2;
        vSSize = vSLength;
        vSArrayNew = vSArray1;
        vFSize = vFLength;
        vFArrayNew = vFArray1;
        for optimisationRound=1:steps2
            maxPositions = [];
            newPositions = [];
            stepSize = ( 1/(steps1-1) ) * 2^(-optimisationRound); 
            for vSIndex = 1:vSSize
                positions = find(estimationPoints2{attackTypeLocalIndex}{vBIndex}(:,1) == vSArrayNew(vSIndex));
                temp = zeros(length(positions),3);
                for d=1:length(positions)
                    temp(d,:) = [EveMSE2{attackTypeLocalIndex}{vBIndex}(positions(d)) , AliceMSE2{attackTypeLocalIndex}{vBIndex}(positions(d)) , estimationPoints2{attackTypeLocalIndex}{vBIndex}(positions(d),2)];
                end
                temp = sortrows(temp,3);
                temp2 = [];
                for d=1:length(temp)
                    if temp(d,1)>EveLimit
                        temp2 = [temp2;temp(d,:)];
                    end
                end
                %if temp2 is empty, this is because final entry gives 0, choose final entry for continuation of algorithm
                if isempty(temp2)
                    temp2 = temp(end,:);
                end
                [~,temp2Index] = min(temp2(:,2));
                maxPositions(vSIndex) = temp2(temp2Index,3);

                if temp2(temp2Index,3) == 0
                    newPositions = [newPositions ; vSArrayNew(vSIndex) , temp2(temp2Index,3) + stepSize ];
                elseif temp2(temp2Index,3) == 1
                    newPositions = [newPositions ; vSArrayNew(vSIndex) , temp2(temp2Index,3) - stepSize ];
                else
                    newPositions = [newPositions ; vSArrayNew(vSIndex) , temp2(temp2Index,3) - stepSize ; vSArrayNew(vSIndex) , temp2(temp2Index,3) + stepSize ];
                end

                if vSIndex>1
                    midpoint = ( vSArrayNew(vSIndex) + vSArrayNew(vSIndex-1) ) /2;
                    limits = sort([maxPositions(vSIndex-1);maxPositions(vSIndex)],1);
                    newPositions2 = limits(1):stepSize:limits(2) ;
                    newPositions3 = [];
                    for e=1:length(newPositions2)
                        newPositions3 = [newPositions3; midpoint , newPositions2(e) ];
                    end
                    newPositions = [newPositions; newPositions3];
                end

            end

            [newPositionsCount ,~] = size(newPositions);
            estimationPoints2{attackTypeLocalIndex}{vBIndex} = [estimationPoints2{attackTypeLocalIndex}{vBIndex} ; newPositions];
            AliceFItemp = zeros(newPositionsCount,1);
            AliceMSEtemp = zeros(newPositionsCount,1);
            EveMSEtemp = zeros(newPositionsCount,1);
            roundCountAveragetemp = zeros(newPositionsCount,1);
            undetectionsAveragetemp = zeros(newPositionsCount,1);
            for newPositionIndex=1:newPositionsCount
                [AliceMSEtemp(newPositionIndex) , EveMSEtemp(newPositionIndex) , roundCountAveragetemp(newPositionIndex) , AliceFItemp(newPositionIndex) ,undetectionsAveragetemp(newPositionIndex) ] = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,nMax,trueValueCell{vBIndex},precision,vB,newPositions(newPositionIndex,1),newPositions(newPositionIndex,2),attackTypeLocal) ;
            end

            AliceFI2{attackTypeLocalIndex}{vBIndex} = [ AliceFI2{attackTypeLocalIndex}{vBIndex} ; AliceFItemp] ;
            AliceMSE2{attackTypeLocalIndex}{vBIndex} = [ AliceMSE2{attackTypeLocalIndex}{vBIndex} ; AliceMSEtemp] ;
            EveMSE2{attackTypeLocalIndex}{vBIndex} = [ EveMSE2{attackTypeLocalIndex}{vBIndex} ; EveMSEtemp] ;
            roundCountAverage2{attackTypeLocalIndex}{vBIndex} = [ roundCountAverage2{attackTypeLocalIndex}{vBIndex} ; roundCountAveragetemp] ;
            undetectionsAverage2{attackTypeLocalIndex}{vBIndex} = [ undetectionsAverage2{attackTypeLocalIndex}{vBIndex} ; undetectionsAveragetemp] ;
        end

        temporaryList = [ AliceMSE2{attackTypeLocalIndex}{vBIndex} , EveMSE2{attackTypeLocalIndex}{vBIndex} , roundCountAverage2{attackTypeLocalIndex}{vBIndex} , AliceFI2{attackTypeLocalIndex}{vBIndex} , estimationPoints2{attackTypeLocalIndex}{vBIndex} , undetectionsAverage2{attackTypeLocalIndex}{vBIndex} ] ;

        if vB==1
            temporaryList2 = [];
            for c=1:length(AliceFI2{attackTypeLocalIndex}{vBIndex})
                if temporaryList(c,5) == 0
                    temporaryList2 = [temporaryList2;temporaryList(c,:)];
                end
            end
            temporaryList = temporaryList2;
            clear temporaryList2
        end


        availableIndexes = find(temporaryList(:,2)>EveLimit);
        hybridList = temporaryList(availableIndexes,:);
        [MmaxValue, MmaxIndex] = min(hybridList(:,1));
        outerMaxM{attackTypeLocalIndex}(vBIndex,1) = vBArray(vBIndex);
        for b=1:7
            outerMaxM{attackTypeLocalIndex}(vBIndex,b+1) = hybridList(MmaxIndex,b);
        end

        fileID = fopen([userpath,sprintf('/MSELimitHybrid_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax,attackTypeTotal(attackTypeLocalIndex),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)],'w');
        fopen(fileID);
        fwrite(fileID,outerMaxM{attackTypeLocalIndex},'double');
        fclose(fileID);

        availableIndexes = find(temporaryList(:,5)==0);
        entangledList = temporaryList(availableIndexes,:);
        availableIndexes2 = find(entangledList(:,2)>EveLimit);
        entangledList2 = entangledList(availableIndexes2,:);
        [EmaxValue,EmaxIndex] = min(entangledList2(:,1));

        outerMaxE{attackTypeLocalIndex}(vBIndex,1) = vBArray(vBIndex);
        for b=1:7
            outerMaxE{attackTypeLocalIndex}(vBIndex,b+1) = entangledList2(EmaxIndex,b);
        end
        fileID = fopen([userpath,sprintf('/MSELimitEntangled_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax,attackTypeTotal(attackTypeLocalIndex),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)],'w');
        fopen(fileID);
        fwrite(fileID,outerMaxE{attackTypeLocalIndex},'double');
        fclose(fileID);

        if vB == 1
            temporaryList(:,5) =1;
        end

        availableIndexes = find(temporaryList(:,5)==1);
        separableList = temporaryList(availableIndexes,:);
        availableIndexes2 = find(separableList(:,2)>EveLimit);
        separableList2 = separableList(availableIndexes2,:);
        [SmaxValue,SmaxIndex] = min(separableList2(:,1));

        outerMaxS{attackTypeLocalIndex}(vBIndex,1) = vBArray(vBIndex);
        for b=1:7
            outerMaxS{attackTypeLocalIndex}(vBIndex,b+1) = separableList2(SmaxIndex,b);
        end

        fileID = fopen([userpath,sprintf('/MSELimitSeparable_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax,attackTypeTotal(attackTypeLocalIndex),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)],'w');
        fopen(fileID);
        fwrite(fileID,outerMaxS{attackTypeLocalIndex},'double');
        fclose(fileID);
    end

end
toc
    

end

