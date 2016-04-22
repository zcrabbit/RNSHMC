%% Compare the result on ePDE
nLeap = 6;
nrep = 5;

%% HMC 2D table
HMCresultTable = zeros(5,5);
nLeapPDE = [4 6 8 10 12];
stepsizes = [0.08 0.16 0.24 0.32 0.40];


for i = 1:length(nLeapPDE)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_2dePDE_' num2str(nLeapPDE(i)) 'nLeap_' num2str(stepsizes(j)) 'stepsize_' '*.mat']);
        data = open(['Results/' Files.name]);
        nSamples = length(data.thetaPosterior);
        ESS = CalculateESS(data.thetaPosterior,nSamples-1);
        HMCresultTable(i,j) = min(ESS)/data.Times;
    end
end

%% HMC rerun result
HMCresultTable = zeros(nrep,6);
Files = dir(['Results/' '*' '_HMC_ePDE_' num2str(nLeap) 'nLeap_' num2str(0.16) '*.mat']);
for i = 1:length(Files)
    data = open(['Results/' Files(i).name]);
    nSamples = length(data.thetaPosterior);
    ESS = CalculateESS(data.thetaPosterior,nSamples-1);
    HMCresultTable(i,:) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
end

HMCresultTable = mean(HMCresultTable);

%% RNSHMC rerun result
RNSHMCresultTable = zeros(nrep,6);

Files = dir(['Results/' '*' '_RNSHMC_ePDE_' num2str(nLeap) 'nLeap_' '*.mat']);
for i = 1:length(Files)
    data = open(['Results/' Files(i).name]);
    nSamples = length(data.thetaPosterior);
    ESS = CalculateESS(data.thetaPosterior,nSamples-1);
    RNSHMCresultTable(i,:) = [data.acprat min(ESS) median(ESS) max(ESS) data.Times/nSamples min(ESS)/data.Times];
end

RNSHMCresultTable = mean(RNSHMCresultTable);


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
    