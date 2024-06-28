%%  First published by Sean William Moore 2024-06, © CC 4.0

%This is computationally intensive code that produces the data for figures 6, 7, 8 and 9 of section 4 in Secure quantum-enhanced measurements on a network of sensors by S. W. Moore and J. A. Dunningham

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285

tic
%%  Single Bob
%reset random number stream
rng('shuffle','twister')
color_order = get(gca,'colororder');

userpath = '/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataSecurity/';

precisionRoot = 10;
precision = 2^precisionRoot;
repetitionTrueValue = 1;%128;
repetitionEachValue = 1;%4;

lN1 = 0:100;
nN1 = numel(lN1);
trueValue1 = 2*pi*rand( repetitionTrueValue , 1 );

FList1 = 0:0.01:1;
nF1 = numel(FList1);

Lambda1 = zeros(nF1,nN1);

for iN= 1:nN1
    for iF=1:nF1
        Lambda1(iF,iN) = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,lN1(iN),trueValue1,precision,1,1,FList1(iF),0);
    end
end

%Store raw data
fileID = fopen([userpath,sprintf('/SecurityLambda1ofN_%d_%d_%d_%d.bin' ,nN1,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
fopen(fileID);
fwrite(fileID,Lambda1,'double');
fclose(fileID);

Lambda1rExpected = zeros(1,nF1);
Lambda1mExpected = zeros(1,nF1);
for iF = 1:nF1
    dr = 1/2*FList1(iF)/2;
    dm = 1/4*FList1(iF)/2;

    G1r = dr*(1-dr).^lN1;
    G1m = dm*(1-dm).^lN1;

    Lambda1rExpected(iF) = sum(Lambda1(iF,:).* G1r);
    Lambda1mExpected(iF) = sum(Lambda1(iF,:).* G1m);

end

figure
hold on
plot(FList1,Lambda1mExpected,'-','color',color_order(3,:),'DisplayName','Measure and resend','LineWidth',2)
plot(FList1,Lambda1rExpected,'--','color',color_order(4,:),'DisplayName','Replace','LineWidth',2)
plot(FList1,Lambda1rExpected-Lambda1mExpected,':','color',color_order(6,:),'DisplayName','Difference','LineWidth',2)
legend('location','northwest')
ylim([0 1])
hold off
xlabel('F','FontSize',18)
ylabel('\Lambda_E','FontSize',24)
box on

%Store results
fileID = fopen([userpath,sprintf('/SecurityLambda1r_%d_%d_%d_%d.bin' ,nN1,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
fopen(fileID);
fwrite(fileID,Lambda1rExpected,'double');
fclose(fileID);

fileID = fopen([userpath,sprintf('/SecurityLambda1m_%d_%d_%d_%d.bin' ,nN1,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
fopen(fileID);
fwrite(fileID,Lambda1mExpected,'double');
fclose(fileID);

%%  Multiple Bobs simulations

lB = [2,4,8,16];
lF = [0:0.01:0.1,0.12:0.02,0.2,0.25,0.3,0.35,0.4:.1:1];
lS = [0:0.001:0.2,0.3:0.1:1];
lN = 0:50;
lE = 1-lS;
lM = 1-lF;

nB = numel(lB);
nN = numel(lN);
nF = numel(lF);
nS = numel(lS);

precisionRoot = 10;
precision = 2^precisionRoot;
repetitionTrueValue = 1;%32;
repetitionEachValue = 1;%4;

%creation of parameter values & declaration of MSE arrays
trueValueCell = cell(nB,1);
LambdaSofN = cell(nB,1);
LambdaEofN = cell(nB,1);
for iB=1:nB
    trueValueCell{iB} = 2*pi*rand( repetitionTrueValue , lB(iB) );
    LambdaSofN{iB} = zeros(nF,nN);
    LambdaEofN{iB} = zeros(nF,nN);
end


LambdaE = cell(nB,1);
LambdaS = cell(nB,1);

%Calcualte data
for iB=1:nB
    vB = lB(iB);
    toc
    tic
    disp(['B = ',num2str(vB)])
    for iN= 1:nN
        for iF=1:nF
            LambdaEofN{iB}(iF,iN) = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,lN(iN),trueValueCell{iB},precision,vB,0,lF(iF),0);
            LambdaSofN{iB}(iF,iN) = networksParameterOptimiserMonteCarloDistibutor(repetitionTrueValue,repetitionEachValue,lN(iN),trueValueCell{iB},precision,vB,1,lF(iF),0);
        end
    end
end

%Store raw data
for iB = 1:nB
    vB = lB(iB);
    fileID = fopen([userpath,sprintf('/SecurityLambdaEofN%dB_%d_%d_%d_%d.bin', vB,nN,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
    fopen(fileID);
    fwrite(fileID,LambdaE{iB},'double');
    fclose(fileID);

    fileID = fopen([userpath,sprintf('/SecurityLambdaSofN%dB_%d_%d_%d_%d.bin', vB,nN,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
    fopen(fileID);
    fwrite(fileID,LambdaS{iB},'double');
    fclose(fileID);
end

toc
tic

%%  Multiple Bobs probability distiribution

for iB = 1:nB
    vB=lB(iB);
    for iS = 1:nS
        vS = lS(iS);
        vE = 1-vS;
        for iF = 1:nF
            vF = lF(iF);
%measure and replace with separable
            dSs = vS* ( 1 - ( 1 - .25*vF/2)^vB );
            dSe = vE* ( .25*vF^vB/2 );
            uS = 1 - dSs - dSe;
            if uS <0
                if uS<0.0001
                    error('uS')
                end
                uS = 0;
            end
            SNTS = uS.^(lN-1).*dSs + uS.^lN.*dSe;
            SNTS(1) = SNTS(1) - dSs/uS;
            if sum (SNTS) > 1.00001
                error('SNTS')
            end
            LambdaS{iB}(iS,iF) = sum( SNTS.*LambdaSofN{iB}(iF,:) );
%measure and replace with entangled
            dEs = vS* ( 1 - ( 1 - .5*vF/2)^vB );
            dEe = vE* ( .25*vF^vB/2 );
            uE = 1 - dEs - dEe;
            if uE <0
                if uE<0.0001
                    error('uE')
                end
                uE = 0;
            end
            SNTE = uE.^(lN-1).*dEs + uE.^lN.*dEe;
            SNTE(1) = SNTE(1) - dEs/uE;
            LambdaE{iB}(iS,iF) = sum( SNTE.*LambdaEofN{iB}(iF,:) );
            if sum (SNTE) > 1.00001
                error('SNTE')
            end
        end
    end
end

%Store results
for iB = 1:nB
    vB = lB(iB);
    fileID = fopen([userpath,sprintf('/SecurityLambdaE%dB_%d_%d_%d_%d.bin', vB,nN,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
    fopen(fileID);
    fwrite(fileID,LambdaE{iB},'double');
    fclose(fileID);

    fileID = fopen([userpath,sprintf('/SecurityLambdaS%dB_%d_%d_%d_%d.bin', vB,nN,precisionRoot,repetitionTrueValue,repetitionEachValue)],'w');
    fopen(fileID);
    fwrite(fileID,LambdaS{iB},'double');
    fclose(fileID);
end

[X,Y] = meshgrid(lF,lS);

figure
for a=1:nB
    subplot(2,2,a)
    colormap default
    mesh(X,Y,LambdaE{a},'FaceColor','flat')
    title(['B = ',num2str(lB(a))],'FontSize',18)
    xlabel('F','FontSize',12)
    ylabel('S','FontSize',12)
    zlabel('\Lambda_E','FontSize',18)
end
colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);

figure
for a=1:nB
    subplot(2,2,a)
    colormap default
    mesh(X,Y,LambdaS{a},'FaceColor','flat')
    title(['B = ',num2str(lB(a))],'FontSize',18)
    xlabel('F','FontSize',12)
    ylabel('S','FontSize',12)
    zlabel('\Lambda_E','FontSize',18)
end
colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);


bottom = -0.5;
top = 0.5;

fig = figure;
colormap cool
for a=1:nB
    subplot(2,2,a)
    mesh(X,Y,LambdaE{a}-LambdaS{a},'FaceColor','flat')
    shading interp;
    clim manual
    clim([bottom top]);
    title(['B = ',num2str(lB(a))],'FontSize',18)
    xlabel('F','FontSize',12)
    ylabel('S','FontSize',12)
    zlabel('\Delta\Lambda_E','FontSize',18)
    zlim([-0.5 0.5])
end
colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);

toc
