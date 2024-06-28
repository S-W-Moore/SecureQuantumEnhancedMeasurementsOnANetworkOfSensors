%%  First published by Sean William Moore 2024-06, © CC 4.0

%This is computationally intensive code that produces the data for figure 5 of section 3 in Secure quantum-enhanced measurements on a network of sensors by S. W. Moore and J. A. Dunningham

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285

%   Comment out unwanted scenarios

nList =  unique(ceil(logspace(0,3.3978,100)));
nList2 = 1:1:100;

precisionRoot = 10;
steps1 = 11;
steps2 = 3;
repetitionTrueValue = 64;
repetitionEachValue = 16;

nMax = [100,500,2500];
vBMax = 6;
EveLimit  = 0.5;
attackTypeTotal = 4;

userpath = '/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataProtocol/';


%%  Extract optimised results for 2500 rounds
Hybrid = cell(3,1);
Separable = cell(3,1);
Entangled = cell(3,1);
for a = 1:3
    for b=1:numel(attackTypeTotal)
        fileID = fopen([userpath,sprintf('/MSELimitHybrid_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Hybrid{a,b} = reshape(tempFileExtraction,[vBMax,8]);


        fileID = fopen([userpath,sprintf('/MSELimitSeparable_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Separable{a,b} = reshape(tempFileExtraction,[vBMax,8]);


        fileID = fopen([userpath,sprintf('/MSELimitEntangled_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Entangled{a,b} = reshape(tempFileExtraction,[vBMax,8]);
    end
end

SList = Hybrid{3}(1:vBMax,6);

FList = Hybrid{3}(1:vBMax,7);


%%  Create data on 5 different scenarios involving. 3 protocols, secure and non-secure
BList = 1:6;
precisionRoot = 10;
repetitionTrueValue = 128;
repetitionEachValue = 8;

tic

outputUserPath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Secure/Hybrid';
disp('Start round 1 of 5')

networksMonteCarloInformationGainMain(nList,BList,SList,FList,precisionRoot,repetitionTrueValue,repetitionEachValue,outputUserPath)

toc
tic
outputUserPath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Secure/Separable';
disp('Start round 2 of 5')

SList = zeros(6,1);
FList = Separable{3}(1:vBMax,7);

networksMonteCarloInformationGainMain(nList2,BList,SList,FList,precisionRoot,repetitionTrueValue,repetitionEachValue,outputUserPath)

outputUserPath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Secure/Entangled';
disp('Start round 3 of 5')

SList = ones(6,1);
FList = Entangled{3}(1:vBMax,7);

networksMonteCarloInformationGainMain(nList2,BList,SList,FList,precisionRoot,repetitionTrueValue,repetitionEachValue,outputUserPath)

toc
tic

SList = ones(6,1);
FList = zeros(6,1);


outputUserPath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Insecure/Separable';
disp('Start round 4 of 5')

networksMonteCarloInformationGainMain(nList,BList,SList,FList,precisionRoot,repetitionTrueValue,repetitionEachValue,outputUserPath)
toc
tic

SList = zeros(6,1);
FList = zeros(6,1);


outputUserPath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Insecure/Entangled';
disp('Start round 5 of 5')

networksMonteCarloInformationGainMain(nList,BList,SList,FList,precisionRoot,repetitionTrueValue,repetitionEachValue,outputUserPath)

toc