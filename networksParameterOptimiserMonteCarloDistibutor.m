%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285


function [mseAlice,mseEve,roundCountAverage,AliceFI,undetectionsAverage] = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,nMax,trueValues,precision,vB,vS,vF,attackType)
%networksParameterOptimiserMonteCalroDistibutor builds the arrays for Monte Carlo simulation using vB, vS & vF then distributes them to MonteCarloLooper for calculation

%input:
%   repetitionTrueValue         number of different true values tested
%   repetitionEachValue         number of repetitions for each true value
%   nMax                        maximum rounds
%   trueValues                  parameter values
%   precision                   data points in likelihood function grid approximations
%   vB                          number of Bobs
%   vS                          separable state rate
%   vF                          fidelity checking rate
%   attackType                  scalar or array, 0 Alice only wihtout attack, 1-4 different attacks by Eve on every round
%           1  -   replace with separable state
%           2  -   measure all in random basis send separable state
%           3  -   replace with entangled state
%           4  -   measure as if entangled then replace with entangled state

%output:
%   mseAlice                    Alice's average circular dispersion equivalent to mean square error
%   mseEve                      scalar or array, Eve's average circular dispersion equivalent to mean square error
%   roundCountAverage           scalar or array, average rounds until eve is detected first time and Alice stops protocol
%   AliceFI                     Alice's Fisher information

%further details:


vE = 1-vS;
vM = 1 -vF;

if numel(unique(attackType)) ~= numel(attackType)
    error('networksParameterOptimiserMonteCalroDistibutor attackType values not a unique set')
end

numberOfTypes = numel(attackType) ;

%declaring some output variables
mseEve = zeros(numel(attackType(attackType~=0)),1);
undetectionsAverage = mseEve;
roundCountAverage = mseEve;
%a counter to allow proper assignment to mseEve and roundCountAverage
attackCount = 1;

for attackTypeIndex=1:numberOfTypes
    if attackType(attackTypeIndex) == 0

        %%  Alice

        %   Occurance probability of any set of meausruements being made together with
        %   index-1 parameters involved (from Alice's (Eve's) point of view).
        indvdlStOccrncPrbbltAc = zeros(vB,1);
        for f=1:vB
            indvdlStOccrncPrbbltAc(f) = vE*vM^f*vF^(vB-f);
        end
        indvdlStOccrncPrbbltAc(1) = indvdlStOccrncPrbbltAc(1) + vS*vM;

        FI = 0;
        for a=1:vB
            FI = FI+indvdlStOccrncPrbbltAc(a)*a*nchoosek(vB-1,a-1)/vB;
        end

        AliceFI = FI ;

        if AliceFI < 0.01/nMax
            mseAlice = 1;
        else

            msrmntNmbrVls = [];
            for f=1:vB
                Pnot0 = 1 - (1-indvdlStOccrncPrbbltAc(f))^nMax;
                %Only count highly likely results for efficient use of RAM and computation time.
                if Pnot0 > 0.1/vB 
                    msrmntNmbrVls = [msrmntNmbrVls,f];
                end
            end
            %Always include vB because even if very rare, will always give information on theta if accessed. And 1 because separable states are fairly liekly to create them.
            msrmntNmbrVls = union(1,union(msrmntNmbrVls,vB));

            %   The number of arrays types worth keeping track of
            [~,cntSz] = size(msrmntNmbrVls);

            msrmtnNmbrCntr = zeros(cntSz);
            for f=1:cntSz
                msrmntNmbrCntr(f) = nchoosek(vB,msrmntNmbrVls(f));
                measurementPositionArray{f} = int8( nchoosek(1:vB,msrmntNmbrVls(f)) );
                msrmntRsltsAc{f} = zeros(5,msrmntNmbrCntr(f));
            end

            msrmntRsltsEv = [];
            dt=whos('measurementPositionArray');
            MB=dt.bytes*9.53674e-7;
            if MB>1000
                disp(['Too much RAM, Alice, vB = ',num2str(vB),' vS = ',num2str(vS), ' vF = ',num2str(vF)])
            end

            mseAlice = networksMonteCarloLooper(repetitionTrueValue,repetitionEachValue,nMax,trueValues,precision,vB,vS,vF,0,msrmntNmbrVls,measurementPositionArray,msrmntNmbrCntr,msrmntRsltsAc,msrmntRsltsEv);

        end

        %%  Eve

    elseif attackType(attackTypeIndex) ==1 || attackType(attackTypeIndex) == 2
        %%  Separable attack
        indvdlStOccrncPrbbltEv = vM;

        EveSFI = 0;
        for a=1:1
            EveSFI = EveSFI+indvdlStOccrncPrbbltEv(a)*a*nchoosek(vB-1,a-1)/vB;
        end

        if EveSFI < 0.01/nMax
            mseEve(attackCount) = 1;
            roundCountAverage(attackCount) = nMax;
        else
            %   The number of Bobs performing a measurement that we track
            msrmntNmbrVls = 1;
            %   The number of arrays types worth keeping track of
            cntSz = 1;

            msrmntNmbrCntr = vB;
            measurementPositionArray{1} = int8( nchoosek(1:vB,1) );
            msrmntRsltsEv{1} = zeros(5,vB);
            msrmntRsltsAc = [];

            dt=whos('measurementPositionArray');
            MB=dt.bytes*9.53674e-7;
            if MB>1000
                disp(['Too much RAM, Alice, vB = ',num2str(vB),' vS = ',num2str(vS), ' vF = ',num2str(vF)])
            end

            [~,mseEve(attackCount),roundCountAverage(attackCount),undetectionsAverage(attackCount)] = networksMonteCarloLooper(repetitionTrueValue,repetitionEachValue,nMax,trueValues,precision,vB,vS,vF,attackType(attackTypeIndex),msrmntNmbrVls,measurementPositionArray,msrmntNmbrCntr,msrmntRsltsAc,msrmntRsltsEv);

            attackCount = attackCount+1;

        end

    elseif attackType(attackTypeIndex) == 3 || attackType(attackTypeIndex) == 4
        %%  Entangled attack

        indvdlStOccrncPrbbltEv = zeros(vB,1);
        for f=1:vB
            indvdlStOccrncPrbbltEv(f) = vM^f*vF^(vB-f);
        end

        EveFI = 0;
        for a=1:vB
            EveFI = EveFI+indvdlStOccrncPrbbltEv(a)*a*nchoosek(vB-1,a-1)/vB;
        end

        if EveFI < 0.01/nMax
            mseEve(attackCount) = 1;
            roundCountAverage(attackCount) = nMax;
        else

            msrmntNmbrVls = [];
            for f=1:vB
                Pnot0 = 1 - (1-indvdlStOccrncPrbbltEv(f))^nMax;

                if Pnot0 > 0.1/vB
                    msrmntNmbrVls = [msrmntNmbrVls,f];
                end
            end
            %Always include vB because even if very rare, will always give information on theta if accessed. And 1 because only need vB of them to get information.
            msrmntNmbrVls = union(1,union(msrmntNmbrVls,vB));

            %   The number of arrays types worth keeping track of
            [~,cntSz] = size(msrmntNmbrVls);

            for f=1:cntSz
                msrmntNmbrCntr(f) = nchoosek(vB,msrmntNmbrVls(f));
                measurementPositionArray{f} = int8( nchoosek(1:vB,msrmntNmbrVls(f)) );
                msrmntRsltsEv{f} = zeros(5,msrmntNmbrCntr(f));
            end

            msrmntRsltsAc = [];
            dt=whos('measurementPositionArray');
            MB=dt.bytes*9.53674e-7;
            if MB>1000
                disp(['Too much RAM, Alice, vB = ',num2str(vB),' vS = ',num2str(vS), ' vF = ',num2str(vF)])
            end

            [~,mseEve(attackCount),roundCountAverage(attackCount),undetectionsAverage(attackCount)] = networksMonteCarloLooper(repetitionTrueValue,repetitionEachValue,nMax,trueValues,precision,vB,vS,vF,attackType(attackTypeIndex),msrmntNmbrVls,measurementPositionArray,msrmntNmbrCntr,msrmntRsltsAc,msrmntRsltsEv);
            attackCount = attackCount+1;
        end

    else
        error('monteCarloSim4 attack type')
    end

end

if isempty(AliceFI)
    AliceFI = NaN;
end

end
