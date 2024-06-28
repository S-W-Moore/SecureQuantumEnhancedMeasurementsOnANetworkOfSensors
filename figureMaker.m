%%  First published by Sean William Moore 2024-06, © CC 4.0

%This code produces the plots in Secure quantum-enhanced measurements on a network of sensors by S. W. Moore and J. A. Dunningham from stored data

%Please cite:
%   Sean William Moore and Jacob Andrew Dunningham. Secure quantum-enhanced measurements on a network of sensors. 2024. arXiv: 2406.19285 [quant-ph]. url: https://arxiv.org/abs/2406.19285

clear

nMax = [100,500,2500];
EveLimit = 0.5;
attackTypeTotal = 4;
color_order = get(gca,'colororder');
vBMax = 6;
vBStep = 1;
vBMin = 1;

precisionRoot = 10;
steps1 = 11;
steps2 = 3;
repetitionTrueValue = 64;
repetitionEachValue = 16;

Hybrid = cell(3,1);
Entangled = cell(3,1);
Separable = cell(3,1);

userpath = '/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataProtocol/';

for a = 1:3
    for b=1:numel(attackTypeTotal)
        fileID = fopen([userpath,sprintf('/MSELimitHybrid_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Hybrid{a,b} = reshape(tempFileExtraction,[vBMax,8]);


        fileID = fopen([userpath,sprintf('/MSELimitEntangled_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Entangled{a,b} = reshape(tempFileExtraction,[vBMax,8]);

        fileID = fopen([userpath,sprintf('/MSELimitSeparable_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Separable{a,b} = reshape(tempFileExtraction,[vBMax,8]);
    end
end


SList = Hybrid{3}(1:vBMax,6);
EList = 1-SList;
FList = Hybrid{3}(1:vBMax,7);
MList = 1-FList;

X=1:vBMax;
X2 = 2:vBMax;

gridHorizontal = 1./[100,500,2500];
for a=1:3
    gridDiagonal(a,:) = X./(10^(a+1));
end
gridDiagonal(1,:) = X./100;
gridDiagonal(2,:) = X./500;
gridDiagonal(3,:) = X./2500;

close all

%%  Figure 2

figure
plot(X,Entangled{1}(1:vBMax,2),'-','color',color_order(3,:),'DisplayName','100 Entangled','LineWidth',2)
hold on
plot(X,Entangled{2}(1:vBMax,2),'--','color',color_order(3,:),'DisplayName','500 Entangled','LineWidth',2)
plot(X,Entangled{3}(1:vBMax,2),':','color',color_order(3,:),'DisplayName','2500 Entangled','LineWidth',2)

plot(X,Separable{1}(1:vBMax,2),'-','color',color_order(4,:),'DisplayName','100 Separable','LineWidth',2)
plot(X,Separable{2}(1:vBMax,2),'--','color',color_order(4,:),'DisplayName','500 Separable','LineWidth',2)
plot(X,Separable{3}(1:vBMax,2),':','color',color_order(4,:),'DisplayName','2500 Separable','LineWidth',2)

plot(X,Hybrid{1}(1:vBMax,2),'-','color',color_order(6,:),'DisplayName','100 Hybrid','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,2),'--','color',color_order(6,:),'DisplayName','500 Hybrid','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,2),':','color',color_order(6,:),'DisplayName','2500 Hybrid','LineWidth',2)

box on
hold off
legend('Location','northwest')
yscale('log')
xlabel('B','FontSize',18)
ylabel('\Lambda_A','FontSize',24)
%ylabel('Mean cirular mean square error of likelihood functions (log)')
xticks(X)
grid on

%%  Figure 3

figure

subplot(2,2,2)
plot(X,Entangled{1}(1:vBMax,4),'-','color',color_order(3,:),'DisplayName','100 Entangled','LineWidth',2)
hold on
plot(X,Entangled{2}(1:vBMax,4),'--','color',color_order(3,:),'DisplayName','500 Entangled','LineWidth',2)
plot(X,Entangled{3}(1:vBMax,4),':','color',color_order(3,:),'DisplayName','2500 Entangled','LineWidth',2)

plot(X,Separable{1}(1:vBMax,4),'-','color',color_order(4,:),'DisplayName','100 Separable','LineWidth',2)
plot(X,Separable{2}(1:vBMax,4),'--','color',color_order(4,:),'DisplayName','500 Separable','LineWidth',2)
plot(X,Separable{3}(1:vBMax,4),':','color',color_order(4,:),'DisplayName','2500 Separable','LineWidth',2)

plot(X,Hybrid{1}(1:vBMax,4),'-','color',color_order(6,:),'DisplayName','100 Hybrid','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,4),'--','color',color_order(6,:),'DisplayName','500 Hybrid','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,4),':','color',color_order(6,:),'DisplayName','2500 Hybrid','LineWidth',2)

box on
hold off
xlabel('B','FontSize',12)
ylabel('Average detection round (log)','FontSize',12)
yscale log
ylim([1 1000])
xticks(X)
grid on
title('(b)','FontSize',18)


subplot(2,2,3)
plot(X,Entangled{1}(1:vBMax,7),'-','color',color_order(3,:),'DisplayName','100 Entangled','LineWidth',2)
hold on
plot(X,Entangled{2}(1:vBMax,7),'--','color',color_order(3,:),'DisplayName','500 Entangled','LineWidth',2)
plot(X,Entangled{3}(1:vBMax,7),':','color',color_order(3,:),'DisplayName','2500 Entangled','LineWidth',2)


plot(X,Separable{1}(1:vBMax,7),'-','color',color_order(4,:),'DisplayName','100 Separable','LineWidth',2)
plot(X,Separable{2}(1:vBMax,7),'--','color',color_order(4,:),'DisplayName','500 Separable','LineWidth',2)
plot(X,Separable{3}(1:vBMax,7),':','color',color_order(4,:),'DisplayName','2500 Separable','LineWidth',2)

plot(X,Hybrid{1}(1:vBMax,7),'-','color',color_order(6,:),'DisplayName','100 Hybrid','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,7),'--','color',color_order(6,:),'DisplayName','500 Hybrid','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,7),':','color',color_order(6,:),'DisplayName','2500 Hybrid','LineWidth',2)

box on
hold off
xlabel('B','FontSize',12)
ylabel('F','FontSize',12)
ylim([0 1])
xticks(X)
title('(c)','FontSize',18)


subplot(2,2,4)
hold on

plot(X2,Hybrid{1}(2:vBMax,6),'-','color',color_order(6,:),'DisplayName','100 Hybrid','LineWidth',2)
plot(X2,Hybrid{2}(2:vBMax,6),'--','color',color_order(6,:),'DisplayName','500 Hybrid','LineWidth',2)
plot(X2,Hybrid{3}(2:vBMax,6),':','color',color_order(6,:),'DisplayName','2500 Hybrid','LineWidth',2)

box on
hold off
xlabel('B','FontSize',12)
ylabel('S','FontSize',12)
ylim([0 1])
xlim([2 vBMax])
xticks(X2)
title('(d)','FontSize',18)

subplot(2,2,1)
plot(X,Entangled{1}(1:vBMax,8),'-','color',color_order(3,:),'DisplayName','100 Entangled','LineWidth',2)
hold on
plot(X,Entangled{2}(1:vBMax,8),'--','color',color_order(3,:),'DisplayName','500 Entangled','LineWidth',2)
plot(X,Entangled{3}(1:vBMax,8),':','color',color_order(3,:),'DisplayName','2500 Entangled','LineWidth',2)

plot(X,Separable{1}(1:vBMax,8),'-','color',color_order(4,:),'DisplayName','100 Separable','LineWidth',2)
plot(X,Separable{2}(1:vBMax,8),'--','color',color_order(4,:),'DisplayName','500 Separable','LineWidth',2)
plot(X,Separable{3}(1:vBMax,8),':','color',color_order(4,:),'DisplayName','2500 Separable','LineWidth',2)

plot(X,Hybrid{1}(1:vBMax,8),'-','color',color_order(6,:),'DisplayName','100 Hybrid','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,8),'--','color',color_order(6,:),'DisplayName','500 Hybrid','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,8),':','color',color_order(6,:),'DisplayName','2500 Hybrid','LineWidth',2)

box on
hold off
xlabel('B','FontSize',12)
ylabel('Undetected proportion','FontSize',12)
ylim([0 1])
xticks(X)
title('(a)','FontSize',18)


fig = gcf;
Lgnd = legend('show');
Lgnd.Position(1) = .13;
Lgnd.Position(2) = .655;

%%  Figure 4

figure
hold on

plot(X,Hybrid{1}(1:vBMax,2),'-','color',color_order(6,:),'DisplayName','100 Hybrid','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,2),'--','color',color_order(6,:),'DisplayName','500 Hybrid','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,2),':','color',color_order(6,:),'DisplayName','2500 Hybrid','LineWidth',2)

for a=1:3
    if a==1
        yline(gridHorizontal(a),'-.k','DisplayName','Cramer-Rao bound, entangled states','LineWidth',2)
    else
        yline(gridHorizontal(a),'-.k','HandleVisibility','off','LineWidth',2)
    end
end

for a=1:3
    if a==1
        plot(X,gridDiagonal(a,:),'-.','DisplayName','Cramer-Rao bound, separable states','color',color_order(2,:),'LineWidth',2)
    else
        plot(X,gridDiagonal(a,:),'-.','HandleVisibility','off','color',color_order(2,:),'LineWidth',2)
    end
end

box on
hold off
legend('Location','northwest')
yscale('log')
xlabel('B','FontSize',18)
ylabel('\Lambda_A','FontSize',24)
%ylabel('Mean cirular mean square error of likelihood functions (log)')
xticks(X)

%%  Figure 5


repetitionTrueValue = 128;
repetitionEachValue = 8;

nList = unique(ceil(logspace(0,3.3978,100)));
BList = 1:6;
BList2 = [1,2,4,6];

nNumber = numel(nList);
BNumber = numel(BList);
BNumber2 = numel(BList2);

userpath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Secure/Hybrid';

fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataHybridSecure = reshape(tempFileExtraction,[nNumber,BNumber]);

userpath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Insecure/Entangled';

fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataEntangled = reshape(tempFileExtraction,[nNumber,BNumber]);

userpath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataInformationGain/Insecure/Separable';

fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataSeparable = reshape(tempFileExtraction,[nNumber,BNumber]);

figure
hold on
for a=1:BNumber2
    a2 = BList2(a);
    subplot(2,2,a)
    hold on
    plot(nList,InformationGainDataEntangled(:,a2),'--','color',color_order(3,:),'DisplayName','Entangled unsecured','LineWidth',2)
    plot(nList,InformationGainDataSeparable(:,a2),':','color',color_order(4,:),'DisplayName','Separable unsecured','LineWidth',2)
    plot(nList,InformationGainDataHybridSecure(:,a2),'color',color_order(6,:),'DisplayName','Hybrid secure','LineWidth',2)
    box on
    hold off
    title(['B =  ',num2str(BList(a2))],'FontSize',18)
    xscale('log')
    xlim([1 2500])
    ylim([1/2500 1])
    xticks([1 10 100 1000])
    yticks([1/1000 1/100 1/10 1])
    yscale('log')
    xlabel('N','FontSize',12)
    ylabel('\Lambda_A','FontSize',18)
end
Lgnd = legend('show');
Lgnd.Position(1) = .4;
Lgnd.Position(2) = .45;

%%  Figure 6

lN1 = 0:100;
nN1 = numel(lN1);

FList1 = 0:0.01:1;
nF1 = numel(FList1);

userpath = '/Users/sean/Documents/PhD/PhD_Matlab/cleanFiles/SecureQuantumEnhancedMeasurementsOnANetworkOfSensors/dataSecurity/';
precisionRoot = 10;
repetitionTrueValue = 128;
repetitionEachValue = 4;

fileID = fopen([userpath,sprintf('/SecurityLambda1r_%d_%d_%d_%d.bin' ,nN1,precisionRoot,repetitionTrueValue,repetitionEachValue)]);
Lambda1rExpected = fread(fileID,'double');

fileID = fopen([userpath,sprintf('/SecurityLambda1m_%d_%d_%d_%d.bin' ,nN1,precisionRoot,repetitionTrueValue,repetitionEachValue)]);
Lambda1mExpected = fread(fileID,'double');


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

%%  Figures 7, 8 & 9

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
repetitionTrueValue = 32;
repetitionEachValue = 4;

LambdaE = cell(nB,1);
LambdaS = cell(nB,1);

for iB = 1:nB
    vB = lB(iB);
    fileID = fopen([userpath,sprintf('/SecurityLambdaE%dB_%d_%d_%d_%d.bin', vB,nN,precisionRoot,repetitionTrueValue,repetitionEachValue)]);
    LambdaE{iB} = fread(fileID,'double');
    LambdaE{iB} = reshape(LambdaE{iB},[nS,nF]);

    fileID = fopen([userpath,sprintf('/SecurityLambdaS%dB_%d_%d_%d_%d.bin', vB,nN,precisionRoot,repetitionTrueValue,repetitionEachValue)]);
    LambdaS{iB} = fread(fileID,'double');
    LambdaS{iB} = reshape(LambdaS{iB},[nS,nF]);
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




%{
%%  Set2
clear Hybrid Entangled Separable


userpath = '/Users/sean/Documents/PhD/PhD_Matlab/Paper02code/data04/2/';

nMax = 100;
EveLimit = 0.5;
attackTypeTotal = [1,2,3,4];%[1,2,3,4]; %at least one of 1 to 4


Hybrid = cell(1,4);
Entangled = cell(1,4);
Separable = cell(1,4);

for a = 1
    for b=1:numel(attackTypeTotal)
        fileID = fopen([userpath,sprintf('/MSELimitHybrid_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Hybrid{a,b} = reshape(tempFileExtraction,[vBMax,8]);


        fileID = fopen([userpath,sprintf('/MSELimitEntangled_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Entangled{a,b} = reshape(tempFileExtraction,[vBMax,8]);

        fileID = fopen([userpath,sprintf('/MSELimitSeparable_L%dnM%dAT%dP%dS%dS%dR%dR%d.bin',EveLimit,nMax(a),attackTypeTotal(b),precisionRoot,steps1,steps2,repetitionTrueValue,repetitionEachValue)]);
        tempFileExtraction = fread(fileID,'double');
        Separable{a,b} = reshape(tempFileExtraction,[vBMax,8]);
    end
end



figure
hold on

plot(X,Hybrid{1}(1:vBMax,2),'-','color',color_order(3,:),'DisplayName','Replace separable','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,2),'--','color',color_order(4,:),'DisplayName','Measure and replace separable','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,2),'-.','color',color_order(6,:),'DisplayName','Replace entangled','LineWidth',2)
plot(X,Hybrid{4}(1:vBMax,2),':','color',color_order(7,:),'DisplayName','Measure and replace entangled','LineWidth',2)


hold off
legend('Location','northwest')
yscale('log')
ylim([10^-3 1])
xlabel('Number of Bobs')
ylabel('Average dispersion of likelihood functions around \theta (log)')
xticks(X)
grid on
title('100 rounds, Hybrid state protocol varying attack types')

figure
hold on 
plot(X,Hybrid{1}(1:vBMax,4),'-','color',color_order(3,:),'DisplayName','Replace separable','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,4),'--','color',color_order(4,:),'DisplayName','Measure and replace separable','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,4),'-.','color',color_order(6,:),'DisplayName','Replace entangled','LineWidth',2)
plot(X,Hybrid{4}(1:vBMax,4),':','color',color_order(7,:),'DisplayName','Measure and replace entangled','LineWidth',2)

hold off
legend('Location','northwest')
xlabel('Number of Bobs')
ylabel('Average rounds before detection (log)')
yscale log
ylim([1 1000])
xticks(X)
grid on
title('100 rounds, Hybrid state protocol varying attack types')


figure
hold on 
plot(X,Hybrid{1}(1:vBMax,7),'-','color',color_order(3,:),'DisplayName','Replace separable','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,7),'--','color',color_order(4,:),'DisplayName','Measure and replace separable','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,7),'-.','color',color_order(6,:),'DisplayName','Replace entangled','LineWidth',2)
plot(X,Hybrid{4}(1:vBMax,7),':','color',color_order(7,:),'DisplayName','Measure and replace entangled','LineWidth',2)

hold off
legend('Location','northwest')
xlabel('Bob number')
ylabel('P_F')
ylim([0 1])
xticks(X)
title('100 rounds, Hybrid state protocol varying attack types')

figure
hold on

plot(X,Hybrid{1}(1:vBMax,6),'-','color',color_order(3,:),'DisplayName','Replace separable','LineWidth',2)
plot(X,Hybrid{2}(1:vBMax,6),'--','color',color_order(4,:),'DisplayName','Measure and replace separable','LineWidth',2)
plot(X,Hybrid{3}(1:vBMax,6),'-.','color',color_order(6,:),'DisplayName','Replace entangled','LineWidth',2)
plot(X,Hybrid{4}(1:vBMax,6),':','color',color_order(7,:),'DisplayName','Measure and replace entangled','LineWidth',2)

hold off
legend('Location','northwest')
xlabel('Bob number')
ylabel('P_S')
ylim([0 1])
xlim([2 vBMax])
xticks(X2)
title('1 000 rounds, Hybrid state protocol varying attack types')
%}

%%  Set3

%{
clear Hybrid Entangled Separable


userpath = '/Users/sean/Documents/PhD/PhD_Matlab/Paper02code/dataInformationGain/Secure/Hybrid';

repetitionTrueValue = 128;%32;%100;%32;
repetitionEachValue = 8;%10;%8;

nList = unique(ceil(logspace(0,3.3978,100)));%unique(floor(logspace(0,4,100)));%[1:1:10,20:10,100,200:100:1000,2000:1000:10000];
BList = 1:6;
BList2 = [1,2,4,6];%1:6;

nNumber = numel(nList);
BNumber = numel(BList);
BNumber2 = numel(BList2);

fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataHybridSecure = reshape(tempFileExtraction,[nNumber,BNumber]);

figure
hold on
for a=1:BNumber
    plot(nList,InformationGainDataHybridSecure(:,a),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
xscale('log')
ylim([10^-4 1])
yscale('log')

title('Hybrid initial states, secure','FontSize',24)


userpath = '/Users/sean/Documents/PhD/PhD_Matlab/Paper02code/dataInformationGain/Insecure/Entangled';



fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataEntangled = reshape(tempFileExtraction,[nNumber,BNumber]);

figure
hold on
for a=1:BNumber
    plot(nList,InformationGainDataEntangled(:,a),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
xscale('log')
ylim([10^-4 1])
yscale('log')
title('Entangled initial states, unsecured','FontSize',24)


userpath = '/Users/sean/Documents/PhD/PhD_Matlab/Paper02code/dataInformationGain/Insecure/Seprabale';



fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataSeparable = reshape(tempFileExtraction,[nNumber,BNumber]);

figure
hold on
for a=1:BNumber2
    a2 = BList2(a);
    plot(nList,InformationGainDataSeparable(:,a),'DisplayName',['B =  ',num2str(BList(a2))])
end
box on
hold off
legend
xscale('log')
ylim([10^-4 1])
yscale('log')
title('Separable initial states, unsecured','FontSize',24)

figure
subplotDimension1 = ceil(sqrt(BNumber));
if subplotDimension1*(subplotDimension1-1) <= BNumber
    subplotDimension2 = subplotDimension1-1;
else
    subplotDimension2 = subplotDimension1;
end
hold on
for a=1:BNumber2
    a2 = BList2(a);
    subplot(2,2,a)
    hold on
    plot(nList,InformationGainDataEntangled(:,a2),'--','color',color_order(3,:),'DisplayName','Entangled unsecured','LineWidth',2)
    plot(nList,InformationGainDataSeparable(:,a2),':','color',color_order(4,:),'DisplayName','Separable unsecured','LineWidth',2)
    plot(nList,InformationGainDataHybridSecure(:,a2),'color',color_order(6,:),'DisplayName','Hybrid secure','LineWidth',2)
    box on
    hold off
    title(['B =  ',num2str(BList(a2))],'FontSize',18)
    %legend
    xscale('log')
    xlim([1 2500])
    ylim([1/2500 1])
    xticks([1 10 100 1000])
    yticks([1/1000 1/100 1/10 1])
    yscale('log')
    xlabel('N','FontSize',12)%number of rounds')
    ylabel('{\Lambda}','FontSize',18)
end
%{
% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
fig.Position(2) = fig.Position(2) + 500;
Lgnd.Position(1) = 0.4;
Lgnd.FontSize = 15;
Lgnd.Position(2) = 0.0;
%}
Lgnd = legend('show');
Lgnd.Position(1) = .4;
Lgnd.Position(2) = .45;


nList2 = 1:1:100;
nNumber2 = 100;

userpath = '/Users/sean/Documents/PhD/PhD_Matlab/Paper02code/data04/3/EvS';
fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataEveSeparable = reshape(tempFileExtraction,[nNumber2,BNumber]);

userpath = '/Users/sean/Documents/PhD/PhD_Matlab/Paper02code/data04/3/EvE';
fileID = fopen([userpath,sprintf('/MSEInformationGainR%dR%d.bin',repetitionTrueValue,repetitionEachValue)]);
tempFileExtraction = fread(fileID,'double');
InformationGainDataEveEntangled = reshape(tempFileExtraction,[nNumber2,BNumber]);


figure
subplot(1,2,1)
hold on
for a=1:BNumber
    plot(nList2,InformationGainDataEveSeparable(:,a),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
%xscale('log')
ylim([10^-2 1])
%yscale('log')
title('Separable initial state','FontSize',24)

subplot(1,2,2)
hold on
for a=1:BNumber
    plot(nList2,InformationGainDataEveEntangled(:,a),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
%xscale('log')
ylim([10^-2 1])
%yscale('log')
title('Entangled initial state','FontSize',24)

nList3 = [0,nList2];
d_r = 1/2;
d_mr = 1/4;

for B=1:6
    d_Sr = SList(B)*(1-(1-d_r*FList(B))^B);
    d_Ssmr = SList(B)*(1-(1-d_mr*FList(B))^B);

    d_Er = EList(B)*d_r*FList(B)^B;
    d_Esmr = EList(B)*d_mr*FList(B)^B;

    u_r = 1 - d_Sr - d_Er;
    u_smr = 1-d_Ssmr - d_Esmr;
    SNT1r(B,:) = u_r.^(nList3-1)*d_Sr + u_r.^(nList3)*d_Er;
    SNT1r(B,1) = SNT1r(B,1) - u_r.^(nList3(1)-1)*d_Sr;
    SNT1smr(B,:) = u_smr.^(nList3-1)*d_Ssmr + u_smr.^(nList3)*d_Esmr;
    SNT1smr(B,1) = SNT1smr(B,1) - u_smr.^(nList3(1)-1)*d_Ssmr ;
end

figure
subplot(1,2,1)
hold on
for a=1:BNumber
    plot(nList3,SNT1r(a,:),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
%xscale('log')
%ylim([10^-2 1])
%yscale('log')
title('Entangled replace','FontSize',24)

subplot(1,2,2)
hold on
for a=1:BNumber
    plot(nList3,SNT1smr(a,:),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
%xscale('log')
%ylim([10^-2 1])
%yscale('log')
title('Separable measure and replace','FontSize',24)

%InformationGainEveEntangled2 = [ones(6,1),InformationGainEveEntangled];
%InformationGainEveSeparable2 = [ones(6,1),InformationGainSeparable];

MSE_r = zeros(BNumber,100);
MSE_smr = zeros(BNumber,100);

for B=1:BNumber
    for measurementNumber=1:100
        if measurementNumber ==1
            MSE_r(B,1) = InformationGainDataEveEntangled(measurementNumber,B);
            MSE_smr(B,1) = InformationGainDataEveSeparable(measurementNumber,B);
        else
            remainingProbability1 = 1 - sum(SNT1r(B,1:measurementNumber)) ;
            remainingProbability2 = 1 - sum(SNT1smr(B,1:measurementNumber)) ;
            MSE_r(B,measurementNumber) = sum( SNT1r(B,2:measurementNumber)*InformationGainDataEveEntangled(1:measurementNumber-1,B) ) + remainingProbability1*InformationGainDataEveEntangled(measurementNumber,B);
            MSE_smr(B,measurementNumber) = sum( SNT1smr(B,2:measurementNumber)*InformationGainDataEveSeparable(1:measurementNumber-1,B) ) + remainingProbability2*InformationGainDataEveSeparable(measurementNumber,B);
        end
    end
end

figure
subplot(1,2,1)
hold on
for a=1:BNumber
    plot(nList2,MSE_r(a,:),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
%xscale('log')
%ylim([10^-2 1])
yscale('log')
title('Entangled replace','FontSize',24)

subplot(1,2,2)
hold on
for a=1:BNumber
    plot(nList2,MSE_smr(a,:),'DisplayName',['B =  ',num2str(BList(a))])
end
box on
hold off
legend
%xscale('log')
%ylim([10^-2 1])
yscale('log')
title('Separable measure and replace','FontSize',24)
%}