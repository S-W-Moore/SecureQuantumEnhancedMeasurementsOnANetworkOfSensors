%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [mseAlice,mseEve,roundCountAverage,averageUndetections] = networksMonteCarloLooper(repetitionTrueValue,repetitionEachValue,nMax,trueValues,precision,vB,vS,vF,attackType,msrmntNmbrVls,measurementPositionArray,msrmntNmbrCntr,msrmntRsltsAc,msrmntRsltsEv)
%networksMonteCarloLooper performs many Monte Caro simulations of the protocol and calls likelihoodOptimiser to get an optimised likelihood function

%input:
%   repetitionTrueValue         number of different true values tested
%   repetitionEachValue         number of repetitions for each true value
%   nMax                        maximum rounds
%   trueValues                  parameter values
%   precision                   data points in likelihood function grid approximations
%   vB                          number of Bobs
%   vS                          separable state rate
%   vF                          fidelity checking rate
%   attackType                  0 Alice only wihtout attack, 1-4 different attacks by Eve on every round
%           1  -   replace with separable state
%           2  -   measure all in random basis send separable state
%           3  -   replace with entangled state
%           4  -   measure as if entangled then replace with entangled state
%   msrmntNmbrVls               the number of results that are likely to occur given vB, vS & vF
%   measurementPositionArray    array for storing above ETC.
%   msrmntNmbrCntr
%   msrmntRslts                 pre-built array for storing measurement results


%output:
%   mseAlice                    Alice's average circular dispersion equivalent to mean square error
%   mseEve                      scalar or array, Eve's average circular dispersion equivalent to mean square error
%   roundCountAverage           scalar or array, average rounds until eve is detected first time and Alice stops protocol

%further details:
%   to be broken down into two functions later

if attackType~=0 && attackType~= 1 && attackType~= 2 && attackType~= 3 && attackType~= 4
    error('attack type invalid')
else
    type = attackType;
end

if attackType == 0
    vR = 0;
else
    vR=1;
end

rng('shuffle','twister')
msrmntRsltsAc2 = msrmntRsltsAc;
msrmntNmbrCntr2 = msrmntNmbrCntr;
clear msrmntRsltsAc msrmntNmbrCntr

vM = 1- vF;
vE = 1 - vS;
[~,cntSz] = size(msrmntNmbrVls);

dPhi = 2*pi/precision;
Phi = linspace(dPhi, 2*pi, precision);

noEavesdropperDetectedOuter = 0;


loopCountOuter = zeros(repetitionTrueValue,repetitionEachValue);
mseOuterA = zeros(repetitionTrueValue,repetitionEachValue);
mseOuterE = zeros(repetitionTrueValue,repetitionEachValue);
parfor repetitionTrueValueCount=1:repetitionTrueValue
    %for repetitionTrueValueCount=1:repetitionTrueValue
    msrmntNmbrCntr = msrmntNmbrCntr2;
    phi = trueValues(repetitionTrueValueCount,:);
    theta = mod(sum(phi),2*pi);
    llhTotalSum = zeros(1,precision);
    llhETotalSum = zeros(1,precision);
    for repetitionOfTrueValueCount=1:repetitionEachValue
        %% Monte-Carlo simulation

        msrmntRsltsEv2 = msrmntRsltsEv;
        msrmntRsltsAc = msrmntRsltsAc2;

        uncategorisedAc = zeros(1,vB);% {};
        uncategorisedEv = zeros(1,vB);% {};

        %
        stop = 0;
        noEavesdropperDetected = 0;%becomes 1 if no eavesdropppper detected

        loopCount = 0;
        while stop == 0

            %%  1

            if rand<vS
                initialState = int8(zeros(1,vB));
                initialSeparable = uint8(1);
                initialEntangled = uint8(0);
                for vBIndex=1:vB
                    initialState(vBIndex) = uint8( randi([0 3]) );
                end
                initialStateTotal = mod(sum(initialState),4); %mod in place of wrapTo2Pi due to not including 2Pi.
            else
                initialSeparable = uint8( 0);
                initialEntangled = uint8( 1);
                initialStateTotal = int8( randi([0 3]) );
            end

            %%  2

            if rand <vR
                attacked = uint8( 1);
            else
                attacked = uint8( 0);
                if initialSeparable == 1
                    arrivalState = initialState;
                end
                arrivalStateTotal = initialStateTotal;
                arrivalEntangled = initialEntangled;
                arrivalSeparable = initialSeparable;
            end

            if attacked == 1 && type ==1
                % replace with separble attack

                arrivalEntangled = 0;
                arrivalSeparable = 1;

                arrivalState = randi([0,3],vB,1);

            elseif attacked == 1 && type == 2
                % measure & resend with separable attack

                arrivalEntangled = 0;
                arrivalSeparable = 1;

                if initialEntangled == 1
                    remainingState = initialStateTotal;
                end
                for vBIndex=1:vB
                    basis = randi([0 1]) ;
                    if initialEntangled == 1 && vBIndex<vB
                        if rand<.5
                            initialState(vBIndex) = basis;
                            arrivalState(vBIndex) = basis;
                            remainingState = mod(remainingState - basis,4);
                        else
                            initialState(vBIndex) = basis + 2;
                            arrivalState(vBIndex) = basis + 2;
                            remainingState = mod(remainingState - basis - 2,4);
                        end
                    elseif initialEntangled == 1 && vBIndex == vB
                        initialState(vBIndex) = remainingState;
                        if mod(remainingState,2) == basis
                            arrivalState(vBIndex) = initialState(vBIndex);
                        else
                            if rand < .5
                                arrivalState(vBIndex) = basis;
                            else
                                arrivalState(vBIndex) = basis + 2;
                            end
                        end


                    elseif initialSeparable == 1
                        if mod(initialState(vBIndex),2) == basis
                            arrivalState(vBIndex) = initialState(vBIndex);
                        else
                            if rand < .5
                                arrivalState(vBIndex) = basis;
                            else
                                arrivalState(vBIndex) = basis + 2;
                            end
                        end
                    else
                        error('measure and replace error')
                    end
                end
                arrivalStateTotal = mod(sum(arrivalState),4);


            elseif attacked == 1 && type == 3
                %   replace with entangled attack

                arrivalEntangled = 1;
                arrivalSeparable = 0;
                arrivalStateTotal = randi([0,3]);
            elseif attacked == 1 && type == 4
                %   measure and replace with entangled attack

                basis = int8( randi([0 1]) );
                arrivalEntangled = 1;
                arrivalSeparable = 0;
                if initialStateTotal == basis || initialStateTotal == basis + 2
                    arrivalStateTotal = initialStateTotal;
                elseif rand<.5
                    arrivalStateTotal = basis;
                else
                    arrivalStateTotal = basis +2;
                end

            elseif attacked == 1
                error('attack type error')
            end

            %%  3

            interactionCount = uint8( zeros(1,vB) );
            phiSet = zeros(1,vB);
            phiSetTotal = 0;
            epsilon = uint8( zeros(1,vB) );
            epsilonTotal = int8( 0 );
            fCheck = int8( zeros(1,vB) );
            mCheck = int8( zeros(1,vB) );
            measurementPositionCounter = int8(zeros(1,vB));
            for h=1:vB
                if rand< vM
                    mCheck(h) = 1;
                    interactionCount(h) = 1;
                    phiSet(h) = phi(h);
                    phiSetTotal = mod(phiSetTotal + phi(h) , 2*pi);
                    epsilon(h) = 0;
                    measurementPositionCounter(h) = h;
                elseif rand <.5
                    fCheck(h) = 1;
                    epsilon(h) = 1;
                    epsilonTotal = mod(epsilonTotal + 1,4);
                else
                    fCheck(h) = 1;
                end
            end
            measurementCount = sum(measurementPositionCounter ~= 0);
            measurementPositionCounter = measurementPositionCounter(measurementPositionCounter~=0);
            fCheckSum = sum(fCheck==1);
            if fCheckSum + measurementCount ~= vB
                error('3')
            end
            if ismemberFastElementwise(measurementCount,msrmntNmbrVls)
                measurementNumberCounterPosition = find(msrmntNmbrVls == measurementCount );
                measurementPosition = ismemberIndexFastRows(  measurementPositionCounter    ,  measurementPositionArray{ measurementNumberCounterPosition } );%[~,measurementPosition] =ismember(  measurementPositionCounter    ,  measurementPositionArray{ measurementNumberCounterPosition } , 'rows' );
                tempUncategorisedAc = 0;
                tempUncategorisedEv = 0;
            elseif measurementCount == 0
                tempUncategorisedAc = 0;
                tempUncategorisedEv = 0;
            elseif measurementCount~=0
                if initialEntangled == 1
                    tempUncategorisedAc = 1;
                    uncategorisedAc(measurementCount) = uncategorisedAc(measurementCount) + 1;
                end
                if arrivalEntangled == 1 && attacked == 1
                    tempUncategorisedEv = 1;
                    uncategorisedEv(measurementCount) = uncategorisedEv(measurementCount) + 1;
                end
            end

            %%  4

            resultArray = int8( zeros(1,vB) );
            resultTotal = int8( 1 );
            if arrivalSeparable == 1
                for i=1:vB
                    argument = ( double(arrivalState(i)) - double(epsilon(i)) ) * pi/2 + phiSet(i);
                    result = int8( rand<.5*(1+cos(argument)) );
                    if result == 0
                        result = -1;
                    end
                    resultArray(i) = result;
                    resultTotal = resultTotal*result;
                end
            else %arrivalEntangled
                finalRand = rand;
                for i=1:vB
                    if i<vB
                        result = int8( 2*randi(2) - 3 );
                    else
                        argumentTotal = mod(( double( arrivalStateTotal ) - double(epsilonTotal)  )*pi/2 + phiSetTotal,2*pi);
                        result = -1+2*int8( finalRand < .5*(1+double(resultTotal)*cos(argumentTotal)));
                    end
                    resultArray(i) = result;
                    resultTotal = resultTotal*result;
                end
            end

            %%  5

            %Alice tests the fidelity checks
            if initialEntangled == 1 && fCheckSum ==vB
                if initialStateTotal == epsilonTotal && resultTotal == -1
                    stop =1;
                elseif initialStateTotal == mod(epsilonTotal+2, 4) && resultTotal == 1
                    stop = 1;
                end
            elseif initialSeparable ==1
                separableFailure = 0;
                for i=1:vB
                    if fCheck(i) == 1
                        if initialState(i) == epsilon(i) && resultArray(i) == -1
                            separableFailure = separableFailure + 1;
                            stop =1;
                        elseif initialState(i) == mod(epsilon(i)+2, 4) && resultArray(i) == 1
                            stop =1;
                            separableFailure = separableFailure + 1;
                        end
                    end
                end
            end

            %%  6

            loopCount = loopCount + 1;
            if loopCount >= nMax %Automatic due to no zero point for vREstimates
                stop =1;
                noEavesdropperDetected = 1;
            end


            %%  7a   Alice interprets results as +-1 & records them
            if attacked == 0
                if initialSeparable == 1
                    results = zeros(4,vB);
                    if separableFailure == 0

                        for vBIndex=1:vB
                            if mCheck(vBIndex) == 1
                                if ( initialState(vBIndex) == 0 && resultArray(vBIndex) == 1 ) || ( initialState(vBIndex) == 2 && resultArray(vBIndex) == -1 )
                                    results(1,vBIndex) = 1;
                                    msrmntRsltsAc{1}(1,vBIndex) = msrmntRsltsAc{1}(1,vBIndex) + 1;
                                elseif ( initialState(vBIndex) == 1 && resultArray(vBIndex) == 1 ) || ( initialState(vBIndex) == 3 && resultArray(vBIndex) == -1 )
                                    results(2,vBIndex) = 1;
                                    msrmntRsltsAc{1}(2,vBIndex) = msrmntRsltsAc{1}(2,vBIndex) + 1;
                                elseif ( initialState(vBIndex) == 0 && resultArray(vBIndex) == -1 ) || ( initialState(vBIndex) == 2 && resultArray(vBIndex) == 1 )
                                    results(3,vBIndex) = 1;
                                    msrmntRsltsAc{1}(3,vBIndex) = msrmntRsltsAc{1}(3,vBIndex) + 1;
                                elseif  ( initialState(vBIndex) == 1 && resultArray(vBIndex) == -1 ) || ( initialState(vBIndex) == 3 && resultArray(vBIndex) == 1 )
                                    results(4,vBIndex) = 1;
                                    msrmntRsltsAc{1}(4,vBIndex) = msrmntRsltsAc{1}(4,vBIndex) + 1;
                                else
                                    error('7a1')
                                end
                                msrmntRsltsAc{1}(5,vBIndex) = msrmntRsltsAc{1}(5,vBIndex) + 1;
                            end
                        end
                    end
                elseif initialEntangled == 1 && tempUncategorisedAc == 0 && measurementCount~=0
                    initialStateTotalNew = mod( initialStateTotal - epsilonTotal ,4);
                    msrmntRsltsAc{measurementNumberCounterPosition}(5,measurementPosition ) = msrmntRsltsAc{measurementNumberCounterPosition}(5,measurementPosition) + 1;
                    if ( initialStateTotalNew == 0 && resultTotal == 1 ) || ( initialStateTotalNew == 2 && resultTotal == -1 )
                        msrmntRsltsAc{measurementNumberCounterPosition}(1,measurementPosition ) = msrmntRsltsAc{measurementNumberCounterPosition}(1,measurementPosition) + 1;
                    elseif ( initialStateTotalNew == 1 && resultTotal == 1 ) || ( initialStateTotalNew == 3 && resultTotal == -1 )
                        msrmntRsltsAc{measurementNumberCounterPosition}(2,measurementPosition ) = msrmntRsltsAc{measurementNumberCounterPosition}(2,measurementPosition) + 1;
                    elseif ( initialStateTotalNew == 0 && resultTotal == -1 ) || ( initialStateTotalNew == 2 && resultTotal == 1 )
                        msrmntRsltsAc{measurementNumberCounterPosition}(3,measurementPosition ) = msrmntRsltsAc{measurementNumberCounterPosition}(3,measurementPosition) + 1;
                    elseif ( initialStateTotalNew == 1 && resultTotal == -1 ) || ( initialStateTotalNew == 3 && resultTotal == 1 )
                        msrmntRsltsAc{measurementNumberCounterPosition}(4,measurementPosition ) = msrmntRsltsAc{measurementNumberCounterPosition}(4,measurementPosition) + 1;
                    else
                        error('7a2')
                    end
                elseif initialEntangled == 1 && tempUncategorisedAc == 1
                elseif initialEntangled == 1 && measurementCount == 0
                else
                    error('7a3')
                end

                %%  7b  Eve interprets her result if she attacked

            elseif attacked == 1
                if type == 1 || type ==2
                    for vBIndex=1:vB
                        if mCheck(vBIndex) == 1
                            if ( arrivalState(vBIndex) == 0 && resultArray(vBIndex) == 1 ) || ( arrivalState(vBIndex) == 2 && resultArray(vBIndex) == -1 )
                                msrmntRsltsEv2{1}(1,vBIndex) = msrmntRsltsEv2{1}(1,vBIndex) + 1;
                            elseif ( arrivalState(vBIndex) == 1 && resultArray(vBIndex) == 1 ) || ( arrivalState(vBIndex) == 3 && resultArray(vBIndex) == -1 )
                                msrmntRsltsEv2{1}(2,vBIndex) = msrmntRsltsEv2{1}(2,vBIndex) + 1;
                            elseif ( arrivalState(vBIndex) == 0 && resultArray(vBIndex) == -1 ) || ( arrivalState(vBIndex) == 2 && resultArray(vBIndex) == 1 )
                                msrmntRsltsEv2{1}(3,vBIndex) = msrmntRsltsEv2{1}(3,vBIndex) + 1;
                            elseif  ( arrivalState(vBIndex) == 1 && resultArray(vBIndex) == -1 ) || ( arrivalState(vBIndex) == 3 && resultArray(vBIndex) == 1 )
                                msrmntRsltsEv2{1}(4,vBIndex) = msrmntRsltsEv2{1}(4,vBIndex) + 1;
                            else
                                error('7b1')
                            end
                            msrmntRsltsEv2{1}(5,vBIndex) = msrmntRsltsEv2{1}(5,vBIndex) + 1;
                        end
                    end
                elseif (type == 3 || type == 4) && tempUncategorisedEv == 0 && measurementCount ~= 0
                    arrivalStateTotalNew = mod( arrivalStateTotal - epsilonTotal ,4);
                    msrmntRsltsEv2{measurementNumberCounterPosition}(5,measurementPosition ) = msrmntRsltsEv2{measurementNumberCounterPosition}(5,measurementPosition) + 1;
                    if ( arrivalStateTotalNew == 0 && resultTotal == 1 ) || ( arrivalStateTotalNew == 2 && resultTotal == -1 )
                        msrmntRsltsEv2{measurementNumberCounterPosition}(1,measurementPosition ) = msrmntRsltsEv2{measurementNumberCounterPosition}(1,measurementPosition) + 1;
                    elseif ( arrivalStateTotalNew == 1 && resultTotal == 1 ) || ( arrivalStateTotalNew == 3 && resultTotal == -1 )
                        msrmntRsltsEv2{measurementNumberCounterPosition}(2,measurementPosition ) = msrmntRsltsEv2{measurementNumberCounterPosition}(2,measurementPosition) + 1;
                    elseif ( arrivalStateTotalNew == 0 && resultTotal == -1 ) || ( arrivalStateTotalNew == 2 && resultTotal == 1 )
                        msrmntRsltsEv2{measurementNumberCounterPosition}(3,measurementPosition ) = msrmntRsltsEv2{measurementNumberCounterPosition}(3,measurementPosition) + 1;
                    elseif ( arrivalStateTotalNew == 1 && resultTotal == -1 ) || ( arrivalStateTotalNew == 3 && resultTotal == 1 )
                        msrmntRsltsEv2{measurementNumberCounterPosition}(4,measurementPosition ) = msrmntRsltsEv2{measurementNumberCounterPosition}(4,measurementPosition) + 1;
                    else
                        error('7b2')
                    end
                elseif (type == 3 || type == 4) && tempUncategorisedEv == 1
                elseif (type == 3 || type == 4) && measurementCount == 0
                else
                    error('7b3')
                end
            else
                error('attacked')
            end

            for c=1:vB
                if vR~= 1 && uncategorisedAc(c) > nchoosek(vB,c);%vB% 10 * nchoosek(vB,c)
                    error('Too many uncategorised Alice')
                end
                if vR ~= 0 && uncategorisedEv(c) > nchoosek(vB,c);%vB%10 * nchoosek(vB,c)
                    error('Too many uncategorised Eve')
                end
            end
        end

        %%  7c

        %Alice
        if vR ~=1
            llh = zeros(cntSz,precision);

            for countIndex=1:cntSz

                if nnz(msrmntRsltsAc{countIndex}(5,:))==msrmntNmbrCntr(countIndex)
                    llh(countIndex,:) = likelihoodOptimiser(msrmntNmbrVls(countIndex),vB,dPhi,Phi,precision,msrmntRsltsAc{countIndex},10);
                elseif nnz(msrmntRsltsAc{countIndex}(5,:))<msrmntNmbrCntr(countIndex)
                    llh(countIndex,:) = ones(1,precision)/(2*pi);
                else
                    error('7c Alice')
                end
            end

            llhTotal = ones(1,precision);
            for countIndex=1:cntSz
                for b=1:precision
                    llhTotal(b) = llhTotal(b) .* llh(countIndex,b);
                end
                if countIndex>=1
                    llhTotal = llhTotal./(sum(sort(llhTotal))*dPhi);
                end
            end
        else
            llhTotal = ones(1,precision)/dPhi;
        end
        llhTotalSum = llhTotalSum + llhTotal/repetitionEachValue;

        %%Eve
        if vR~=0
            if noEavesdropperDetected == 1
                noEavesdropperDetectedOuter = noEavesdropperDetectedOuter +1;
            end
            llhE = zeros(cntSz,precision);

            for countIndex=1:cntSz
                if nnz(msrmntRsltsEv2{countIndex}(5,:))==msrmntNmbrCntr(countIndex)
                    llhE(countIndex,:) = likelihoodOptimiser(msrmntNmbrVls(countIndex),vB,dPhi,Phi,precision,msrmntRsltsEv2{countIndex},10);
                elseif nnz(msrmntRsltsEv2{countIndex}(5,:)) < msrmntNmbrCntr(countIndex)
                    llhE(countIndex,:) = ones(1,precision)/(2*pi);
                else
                    error('7c Eve')
                end
            end

            llhETotal = ones(1,precision);
            for countIndex=1:cntSz
                for precisionIndex = 1:precision
                    llhETotal(precisionIndex) = llhETotal(precisionIndex) .* llhE(countIndex,precisionIndex);
                end
                if countIndex>=1
                    llhETotal = llhETotal./(sum(sort(llhETotal))*dPhi);
                end
            end
        else
            llhETotal = ones(1,precision);
        end
        llhETotalSum = llhETotalSum + llhETotal/repetitionEachValue;


        [~,~,mseOuterA(repetitionTrueValueCount,repetitionOfTrueValueCount)] = circStats(Phi ,llhTotal , theta);
        [~,~,mseOuterE(repetitionTrueValueCount,repetitionOfTrueValueCount)] = circStats(Phi ,llhETotal , theta);

        loopCountOuter( repetitionTrueValueCount, repetitionOfTrueValueCount) = loopCount;
    end

end

roundCountAverage = mean(mean(loopCountOuter));

mseAlice = mean(mseOuterA,"all");
mseEve = mean(mseOuterE,"all");

averageUndetections = noEavesdropperDetectedOuter / (repetitionEachValue*repetitionTrueValue);
end

function output = ismemberFastElementwise(A,B)

output = ~min(abs( A - B) );

end

function output = ismemberIndexFastRows(A,B)

output = find(abs( sum( A - B , 2 ) )  == 0 );

end