%% compare sampling efficiency with different Leapfrog steps

%% ePDE 2Dsimulation vs stepsize + nLeap
nLeapPDE = [4 6 8 10 12];
stepsizes = [0.08 0.16 0.24 0.32 0.40];
nrep = 5;

for nr = 1:nrep
    for nL = nLeapPDE
        for stepsize = stepsizes
            [thetaPosterior,Times,acprat] = ePDE_HMC(stepsize,nL);
        end
    end
end

%% HMC 2D table
HMCresultTable = zeros(5,5,nrep);
nLeapPDE = [4 6 8 10 12];
stepsizes = [0.08 0.16 0.24 0.32 0.40];


for i = 1:length(nLeapPDE)
    for j = 1:length(stepsizes)
        Files = dir(['Results/' '*' '_HMC_ePDE_' num2str(nLeapPDE(i)) 'nLeap_' num2str(stepsizes(j)) 'stepsize_' '*.mat']);
        for k = 1:length(Files)
            data = open(['Results/' Files(k).name]);
            nSamples = length(data.thetaPosterior);
            ESS = CalculateESS(data.thetaPosterior,nSamples-1);
            HMCresultTable(i,j,k) = min(ESS)/data.Times;
        end
    end
end

HMCresultTable = mean(HMCresultTable,3);
% 2D plot
imagesc(stepsizes,nLeapPDE,HMCresultTable);
set(gca,'Ydir','Normal');