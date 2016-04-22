%% compare sampling efficiency with different Leapfrog steps
%% Simulated 50D
load simulated50d/lrdata
%nLeapSimu = [3 6 12 18 24 30];
% stepsize = 0.045
nLeapSimu = [2 5 10 20 30];
stepsizes = [0.01 0.025 0.05 0.075 0.1];
nrep = 1;

for nr = 1:nrep
    for nL = nLeapSimu
        for stepsize = stepsizes
            [betaPosterior,Times,acprat] = BLR_HMC(y,X,stepsize,nL,'2dSimulation');
        end
    end
end

%% a9a adult
% stepsize = 0.012
load Data/a9a60d
nLeapA9a = [5,10,20,30,50];
stepsizes = [0.004 0.008 0.012 0.016 0.020]; 
nrep = 1;

for nr = 1:nrep
    for nL = nLeapA9a
        for stepsize = stepsizes
            [betaPosterior,Times,acprat] = BLR_HMC(y,X,stepsize,nL,'2da9a60d');
        end
    end
end

%% bank marketing 
load bankmarketing/lrbanknorm
% stepsize = 0.012
nLeapBank = [15 30 45 60 75 90];
stepsizes = [0.004 0.008 0.012 0.014 0.015];
nrep = 5;

for nr = 1:nrep
    for nL = 90
        for stepsize = stepsizes
            [betaPosterior,Times,acprat] = BLR_HMC(y,X,stepsize,nL,'2dBank');
        end
    end
end



%% plots
% 50D simulation
SimuESSTable = zeros(nrep,length(nLeapSimu));
for i = 1:length(nLeapSimu)
    Files = dir(['Results/' '*' 'BLRSimulation_' num2str(nLeapSimu(i)) 'nLeap_' '*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        SimuESSTable(j,i) = min(CalculateESS(data.betaPosterior,length(data.betaPosterior)-1))/data.Times;
    end
end

SimuResults = quantile(SimuESSTable,[0.1 0.5 0.9]);
LeHMCSimu = errorbar(nLeapSimu,SimuResults(2,:),SimuResults(2,:)-SimuResults(1,:),SimuResults(3,:)-SimuResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeHMCSimu,'markerfacecolor',get(LeHMCSimu,'color'));
h = legend(LeHMCSimu,'HMC');
set(h,'FontSize',12);
xlabel('Number of leap-frog steps');
ylabel('Minimum ESS per second');
title('Simulation');
xlim([0,30]);
ylim([0,16]);

% a9a
A9aESSTable = zeros(nrep,length(nLeapA9a));
for i = 1:length(nLeapA9a)
    Files = dir(['Results/' '*' '_HMC_BLRa9a60d_' num2str(nLeapA9a(i)) 'nLeap_' '*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        A9aESSTable(j,i) = min(CalculateESS(data.betaPosterior,length(data.betaPosterior)-1))/data.Times;
    end
end

A9aResults = quantile(A9aESSTable,[0.1 0.5 0.9]);
figure(2);
LeHMCA9a = errorbar(nLeapA9a,A9aResults(2,:),A9aResults(2,:)-A9aResults(1,:),A9aResults(3,:)-A9aResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeHMCA9a,'markerfacecolor',get(LeHMCA9a,'color'));
h = legend(LeHMCA9a,'HMC');
set(h,'FontSize',12);
title('Adult');
xlabel('Number of leap-frog steps');
ylabel('Minimum ESS per second');
xlim([0,40]);


% bank marketing
BankESSTable = zeros(nrep,length(nLeapBank));
for i = 1:length(nLeapA9a)
    Files = dir(['Results/' '*' 'BLRBank_' num2str(nLeapBank(i)) 'nLeap_' '*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        BankESSTable(j,i) = min(CalculateESS(data.betaPosterior,length(data.betaPosterior)-1))/data.Times;
    end
end

BankResults = quantile(BankESSTable,[0.1 0.5 0.9]);
figure(3);
LeHMCBank = errorbar(nLeapBank,BankResults(2,:),BankResults(2,:)-BankResults(1,:),BankResults(3,:)-BankResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeHMCBank,'markerfacecolor',get(LeHMCBank,'color'));
h = legend(LeHMCBank,'HMC');
set(h,'FontSize',12);
title('Bank Marketing');
xlabel('Number of leap-frog steps');
ylabel('Minimum ESS per second');
xlim([20,90]);