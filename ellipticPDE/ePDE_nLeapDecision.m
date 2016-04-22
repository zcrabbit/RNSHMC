%% compare sampling efficiency with different Leapfrog steps

%% ePDE
nLeapPDE = [2 4 6 8 10];
nrep = 5;
D = 20;

for nr = 1:nrep
    for nL = nLeapPDE
        [thetaPosterior,Times,acprat] = ePDE_HMC(D,0.24,nL);
    end
end

%% ePDE 2Dsimulation vs stepsize + nLeap
nLeapPDE = [2 4 6 8 10 12];
stepsizes = [0.08 0.16 0.24 0.32 0.40];
nrep = 1;
D = 20;

for nr = 1:nrep
    for nL = 12
        for stepsize = stepsizes
            [thetaPosterior,Times,acprat] = ePDE_HMC(D,stepsize,nL);
        end
    end
end

%% plots
PDEESSTable = zeros(nrep,length(nLeapPDE));
for i = 1:length(nLeapPDE)
    Files = dir(['Results/' '*' 'ePDE_' num2str(nLeapPDE(i)) 'nLeap_' '*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        PDEESSTable(j,i) = min(CalculateESS(data.thetaPosterior,length(data.thetaPosterior)-1))/data.Times;
    end
end

figure(4);
PDEResults = quantile(PDEESSTable,[0.1 0.5 0.9]);
LeHMCPDE = errorbar(nLeapPDE,PDEResults(2,:),PDEResults(2,:)-PDEResults(1,:),PDEResults(3,:)-PDEResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeHMCPDE,'markerfacecolor',get(LeHMCPDE,'color'));
h = legend(LeHMCPDE,'HMC');
set(h,'FontSize',12);
xlabel('Number of leap-frog steps');
ylabel('Minimum ESS per second');
title('ePDE');
xlim([0,12]);
%ylim([0,16]);