%% Compare the result on ePDE
nLeap = 10;
nrep = 5;


%% HMC
HMCresultTable = zeros(nrep,6);

Files = dir(['Results/' '*' '_HMC_ePDE_' num2str(nLeap) 'nLeap_' '*.mat']);
for i = 1:length(Files)
    data = open(['Results/' Files(i).name]);
    nSamples = length(data.thetaPosterior);
    ESS = CalculateESS(data.thetaPosterior,nSamples-1);
    HMCresultTable(i,:) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
end

HMCresultTable = mean(HMCresultTable);

%% RNSHMC
RNSHMCresultTable = zeros(nrep,6);

Files = dir(['Results/' '*' '_RNSHMC_ePDE_' num2str(nLeap) 'nLeap_' '*.mat']);
for i = 1:length(Files)
    data = open(['Results/' Files(i).name]);
    nSamples = length(data.thetaPosterior);
    ESS = CalculateESS(data.thetaPosterior,nSamples-1);
    RNSHMCresultTable(i,:) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
end

RNSHMCresultTable = mean(RNSHMCresultTable);
    