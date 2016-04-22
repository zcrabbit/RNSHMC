%% compare sampling efficiency with different Leapfrog steps and step sizes
%% Simulated 50D
% stepsize = 0.045
nLeapSimu = [2 5 10 20 30];
stepsizes = [0.01 0.025 0.05 0.075 0.1];
nrep = 5;

for nr = 1:nrep
    for nL = nLeapSimu
        for stepsize = stepsizes
            [betaPosterior,Times,acprat] = BLR_HMC(stepsize,nL,'simulation');
        end
    end
end

%% a9a adult
% stepsize = 0.012
nLeapA9a = [5,10,20,30,50];
stepsizes = [0.004 0.008 0.012 0.016 0.020]; 
nrep = 5;

for nr = 1:nrep
    for nL = nLeapA9a
        for stepsize = stepsizes
            [betaPosterior,Times,acprat] = BLR_HMC(stepsize,nL,'a9a60d');
        end
    end
end

%% bank marketing 
% stepsize = 0.012
nLeapBank = [30 45 60 75 90];
stepsizes = [0.004 0.008 0.012 0.014 0.015];
nrep = 5;

for nr = 1:nrep
    for nL = nLeapBank
        for stepsize = stepsizes
            [betaPosterior,Times,acprat] = BLR_HMC(stepsize,nL,'bankmarket');
        end
    end
end

%% HMC Simulation 2D plot
HMCresultTable = zeros(5,5,nrep);
stepsizes = [0.01 0.025 0.05 0.075 0.1];

for i = 1:length(nLeapSimu)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_BLR' 'simulation' '_' num2str(nLeapSimu(i)) 'nLeap' '_' num2str(stepsizes(j)) 'stepsize' '*.mat']);
        for k = 1:length(Files)
            data = open(['Results/' Files(k).name]);
            nSamples = length(data.betaPosterior);
            ESS = CalculateESS(data.betaPosterior,nSamples-1);
            HMCresultTable(i,j,k) = min(ESS)/data.Times;
        end
    end
end

HMCresultTable = mean(HMCresultTable,3);
% 2D plot
imagesc(stepsizes,nLeapSimu,HMCresultTable);
set(gca,'Ydir','Normal');


%% HMC bank 2D plot
HMCresultTable = zeros(5,5,nrep);
stepsizes = [0.004 0.008 0.012 0.014 0.015];

for i = 1:length(nLeapBank)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_BLR' 'bankmarket' '_' num2str(nLeapBank(i)) 'nLeap' '_' num2str(stepsizes(j)) 'stepsize' '*.mat']);
        for k = 1:length(Files)
            data = open(['Results/' Files(k).name]);
            nSamples = length(data.betaPosterior);
            ESS = CalculateESS(data.betaPosterior,nSamples-1);
            HMCresultTable(i,j,k) = min(ESS)/data.Times;
        end
    end
end

HMCresultTable = mean(HMCresultTable,3);
% 2D plot
imagesc(stepsizes,nLeapBank,HMCresultTable);
set(gca,'Ydir','Normal');


%% HMC a9a 2D plot
HMCresultTable = zeros(5,5,5);
nLeapA9a = [5,10,20,30,50];
stepsizes = [0.004 0.008 0.012 0.016 0.020]; 

for i = 1:length(nLeapA9a)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_BLR' 'a9a60d' '_' num2str(nLeapA9a(i)) 'nLeap' '_' num2str(stepsizes(j)) 'stepsize' '*.mat']);
        for k = 1:length(Files)
            data = open(['Results/' Files(k).name]);
            nSamples = length(data.betaPosterior);
            ESS = CalculateESS(data.betaPosterior,nSamples-1);
            HMCresultTable(i,j) = min(ESS)/data.Times;
        end
    end
end

HMCresultTable = mean(HMCresultTable,3);
% 2D plot
imagesc(stepsizes,nLeapA9a,HMCresultTable);
set(gca,'Ydir','Normal');