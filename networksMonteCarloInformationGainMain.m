%%  First published by Sean William Moore 2024-06, © CC 4.0

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285

function networksMonteCarloInformationGainMain(nList,BList,SList,FList,precisionRoot,repetitionTrueValue,repetitionEachValue,outputUserPath)
%networksMonteCarloInformationGainMain calls networksParameterOptimiserMonteCarloDistibutor repeatedly and stores the values for the n and B inputs

%input:
%   nList                       array of values for number of rounds
%   BList                       array of values for number of Bobs
%   SList                       array of S values corresponding to B values of same index
%   FList                       array of F values corresponding to B values of same index
%   precision                   data points in likelihood function grid approximations is 2^precisionRoot
%   repetitionTrueValue         number of different true values tested
%   repetitionEachValue         number of repetitions for each true value
%   outputUserPath              where to store data

if max(BList)>127
    error('networksMonteCarloInformationGainMain optimisation requires vB to be less than int8')
end

nNumber = numel(nList);
BNumber = numel(BList);

%define output userpath
userpath(outputUserPath)

%simulation parameters to keep
precision = 2^precisionRoot;

%reset random number stream
rng('shuffle','twister')

%creation of parameter values
trueValueCell = cell(BNumber,1);
for BIndx=1:BNumber
    trueValueCell{BIndx} = 2*pi*rand( repetitionTrueValue , BList(BIndx) );
end

%declaration of Round 1 storage variables
AliceMSE = zeros(nNumber,BNumber);

%begin timer for progress tracking
toc

%Call networksParameterOptimiserMonteCarloDistibutor repeatedly and store values
for nIndx=1:nNumber
    toc
    tic
    disp(['round = ',num2str(nIndx) ,' of ',num2str(nNumber),', n = ',num2str(nList(nIndx))])
    for BIndx=1:BNumber
        AliceMSE(nIndx,BIndx) = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,nList(nIndx),trueValueCell{BIndx},precision,BList(BIndx),SList(BIndx),FList(BIndx),0);
    end
    fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)],'w');
    fopen(fileID);
    fwrite(fileID,AliceMSE,'double');
    fclose(fileID);
end
tic

end

