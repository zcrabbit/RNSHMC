%% Compare the results on different data
nLeap = [6,10,45];
Data = {'Simulation' 'a9a60d' 'Bank'};
nrep = 5;






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