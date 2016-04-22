%% Compare the results on different data
nLeap = [6,10,45];
Data = {'Simulation' 'a9a60d' 'Bank'};
nrep = 5;

%% HMC bank 2D table
HMCresultTable = zeros(5,5,5);
nLeapBank = [30 45 60 75 90];
stepsizes = [0.004 0.008 0.012 0.014 0.015];

for i = 1:length(nLeapBank)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_BLR' '2dBank' '_' num2str(nLeapBank(i)) 'nLeap' '_' num2str(stepsizes(j)) 'stepsize' '*.mat']);
        for k = 1:length(Files)
            data = open(['Results/' Files(k).name]);
            nSamples = length(data.betaPosterior);
            ESS = CalculateESS(data.betaPosterior,nSamples-1);
            HMCresultTable(i,j,k) = min(ESS)/data.Times;
        end
    end
end

HMCresultTable = mean(HMCresultTable,3);

%% HMC a9a 2D table
HMCresultTable = zeros(5,5);
nLeapA9a = [5,10,20,30,50];
stepsizes = [0.004 0.008 0.012 0.016 0.020]; 

for i = 1:length(nLeapA9a)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_BLR' '2da9a60d' '_' num2str(nLeapA9a(i)) 'nLeap' '_' num2str(stepsizes(j)) 'stepsize' '*.mat']);     
        data = open(['Results/' Files.name]);
        nSamples = length(data.betaPosterior);
        ESS = CalculateESS(data.betaPosterior,nSamples-1);
        HMCresultTable(i,j) = min(ESS)/data.Times;
    end
end

%% HMC Simulation 2D table
HMCresultTable = zeros(5,5);
nLeapSimu = [2 5 10 20 30];
stepsizes = [0.01 0.025 0.05 0.075 0.1];

for i = 1:length(nLeapSimu)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_BLR' '2dSimulation' '_' num2str(nLeapSimu(i)) 'nLeap' '_' num2str(stepsizes(j)) 'stepsize' '*.mat']);     
        data = open(['Results/' Files.name]);
        nSamples = length(data.betaPosterior);
        ESS = CalculateESS(data.betaPosterior,nSamples-1);
        HMCresultTable(i,j) = min(ESS)/data.Times;
    end
end


%% HMC
HMCresultTable = zeros(length(Data),6,5);

for i = 1:length(Data)
    Files = dir(['Results/' '*' '_HMC_BLR' Data{i} '_' num2str(nLeap(i)) 'nLeap' '*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        nSamples = length(data.betaPosterior);
        ESS = CalculateESS(data.betaPosterior,nSamples-1);
        HMCresultTable(i,:,j) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
    end
end

HMCresultTable = mean(HMCresultTable,3);


%% RNS-HMC
RNSHMCresultTable = zeros(length(Data),6,5);
Data = {'simulation' 'a9a60d' 'bankmarket'};
for i = 1:length(Data)
    Files = dir(['Results/' '*' '_RNSHMC_BLR' Data{i} '_' num2str(nLeap(i)) 'nLeap' '*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        nSamples = length(data.betaPosterior);
        ESS = CalculateESS(data.betaPosterior,nSamples-1);
        RNSHMCresultTable(i,:,j) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
    end
end

RNSHMCresultTable = mean(RNSHMCresultTable,3);

%% RMHMC
nLeap = 4;
Data = 'a9a60d';
RMHMCresultTable = zeros(5,6);

Files = dir(['Results/' '*' '_RMHMC_BLR_' Data '_' num2str(nLeap) 'nLeap_' '*.mat']);
for i = 1:length(Files)
    data = open(['Results/' Files(i).name]);
    nSamples = length(data.betaPosterior);
    ESS = CalculateESS(data.betaPosterior,nSamples-1);
    RMHMCresultTable(i,:) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
end

RMHMCresultTable = mean(RMHMCresultTable);
        

%% RNS-RMHMC
nLeap = 4;
Data = 'a9a60d';
%RNSRMHMCresultTable = zeros(5,6);

Files = dir(['Results/' '*' '_RNSRMHMC_BLR_' Data '_' num2str(nLeap) 'nLeap_' '*.mat']);
for i = 1:length(Files)
    data = open(['Results/' Files(i).name]);
    nSamples = length(data.betaPosterior);
    ESS = CalculateESS(data.betaPosterior,nSamples-1);
    RNSRMHMCresultTable(i,:) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
end

RNSRMHMCresultTable = mean(RNSRMHMCresultTable);