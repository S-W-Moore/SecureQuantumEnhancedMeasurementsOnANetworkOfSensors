%%  First published by Sean William Moore 2024-06, © CC 4.0

%This is computationally intensive code that produces the data for figures 2, 3 and 4 of section 2 in Secure quantum-enhanced measurements on a network of sensors by S. W. Moore and J. A. Dunningham

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285



nMax = [100,500,2500];
EveLimit = 0.5;
attackTypeTotal = 4;

vBMax = 6;
vBStep = 1;
vBMin = 1;

precisionRoot = 10;
steps1 = 11;
steps2 = 3;
repetitionTrueValue = 1;%64;
repetitionEachValue = 1;%16;

outputUserPath = '/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataProtocol/';

for iN = 1:3
    networksMonteCarloOptimisationMain(nMax(iN),EveLimit,attackTypeTotal,vBMin,vBMax,vBStep,precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue,outputUserPath);
end
