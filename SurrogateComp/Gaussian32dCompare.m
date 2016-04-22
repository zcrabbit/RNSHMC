%% Summarize the results for the 32D Gaussian model.

ndata = [300 600 900 1500 2100 3000];
nrep = 10;
%% Results for GP
GPTrainTime = zeros(nrep,length(ndata));
GPAcprat = zeros(nrep,length(ndata));
GPMeanESSperSec = zeros(nrep,length(ndata));

for i = 1:length(ndata)
    Files = dir(['Results/GP_' num2str(ndata(i)) 'n*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        GPTrainTime(j,i) = data.TrainTime;
        GPAcprat(j,i) = data.acprat;
        GPMeanESSperSec(j,i) = mean(CalculateESS(data.thetasaved,length(data.thetasaved)-1))/data.Time;
    end
end
GPMeanESSperSec(isnan(GPMeanESSperSec)) = 0;

%% Results for Additive NN
AddNetTrainTime = zeros(nrep,length(ndata));
AddNetAcprat = zeros(nrep,length(ndata));
AddNetMeanESSperSec = zeros(nrep,length(ndata));

for i = 1:length(ndata)
    Files = dir(['Results/NN_' num2str(ndata(i)) 'n*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        AddNetTrainTime(j,i) = data.TrainTime;
        AddNetAcprat(j,i) = data.acprat;
        AddNetMeanESSperSec(j,i) = mean(CalculateESS(data.thetasaved,length(data.thetasaved)-1))/data.Time;
    end
end


%% Results for FITC
FITCTrainTime = zeros(nrep,length(ndata));
FITCAcprat = zeros(nrep,length(ndata));
FITCMeanESSperSec = zeros(nrep,length(ndata));

for i = 1:length(ndata)
    Files = dir(['Results/FITC_' num2str(ndata(i)) 'n*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        FITCTrainTime(j,i) = data.TrainTime;
        FITCAcprat(j,i) = data.acprat;
        FITCMeanESSperSec(j,i) = mean(CalculateESS(data.thetasaved,length(data.thetasaved)-1))/data.Time;
    end
end

FITCMeanESSperSec(isnan(FITCMeanESSperSec)) = 0;

%% Results for RBF NN
RBFNetTrainTime = zeros(nrep,length(ndata));
RBFNetAcprat = zeros(nrep,length(ndata));
RBFNetMeanESSperSec = zeros(nrep,length(ndata));

for i = 1:length(ndata)
    Files = dir(['Results/RBF_' num2str(ndata(i)) 'n*.mat']);
    for j = 1:length(Files)
        data = open(['Results/' Files(j).name]);
        RBFNetTrainTime(j,i) = data.TrainTime;
        RBFNetAcprat(j,i) = data.acprat;
        RBFNetMeanESSperSec(j,i) = mean(CalculateESS(data.thetasaved,length(data.thetasaved)-1))/data.Time;
    end
end

%% Plots
% TrainTime
figure;
GPResults = quantile(GPTrainTime,[0.10 0.5 0.90]);
LeGP = errorbar(ndata,GPResults(2,:),GPResults(2,:)-GPResults(1,:), GPResults(3,:)-GPResults(2,:),'v--','linewidth',1,'markersize',4);
set(LeGP,'markerfacecolor',get(LeGP,'color'));
hold on;
AddNetResults = quantile(AddNetTrainTime,[0.10 0.5 0.90]);
LeAddNet = errorbar(ndata,AddNetResults(2,:),AddNetResults(2,:)-AddNetResults(1,:), AddNetResults(3,:)-AddNetResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeAddNet,'markerfacecolor',get(LeAddNet,'color'));
FITCResults = quantile(FITCTrainTime,[0.10 0.5 0.90]);
LeFITC = errorbar(ndata,FITCResults(2,:),FITCResults(2,:)-FITCResults(1,:), FITCResults(3,:)-FITCResults(2,:),'o--','linewidth',1,'markersize',4);
set(LeFITC,'markerfacecolor',get(LeFITC,'color'));
RBFNetResults = quantile(RBFNetTrainTime,[0.10 0.5 0.90]);
LeRBFNet = errorbar(ndata,RBFNetResults(2,:),RBFNetResults(2,:)-RBFNetResults(1,:), RBFNetResults(3,:)-RBFNetResults(2,:),'d--','linewidth',1,'markersize',4);
set(LeRBFNet,'markerfacecolor',get(LeRBFNet,'color'));

legend([LeGP LeAddNet LeFITC LeRBFNet],'Full GP','Additive RB','FITC','RBF RB');
xlabel('Number of observations');
ylabel('Training time in seconds');
xlim([0,3000]);


% Acceptance rate
figure;
GPResults = quantile(GPAcprat,[0.10 0.5 0.90]);
LeGP = errorbar(ndata,GPResults(2,:),GPResults(2,:)-GPResults(1,:), GPResults(3,:)-GPResults(2,:),'v--','linewidth',1,'markersize',4);
set(LeGP,'markerfacecolor',get(LeGP,'color'));
hold on;
AddNetResults = quantile(AddNetAcprat,[0.10 0.5 0.90]);
LeAddNet = errorbar(ndata,AddNetResults(2,:),AddNetResults(2,:)-AddNetResults(1,:), AddNetResults(3,:)-AddNetResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeAddNet,'markerfacecolor',get(LeAddNet,'color'));
FITCResults = quantile(FITCAcprat,[0.10 0.5 0.90]);
LeFITC = errorbar(ndata,FITCResults(2,:),FITCResults(2,:)-FITCResults(1,:), FITCResults(3,:)-FITCResults(2,:),'o--','linewidth',1,'markersize',4);
set(LeFITC,'markerfacecolor',get(LeFITC,'color'));
RBFNetResults = quantile(RBFNetAcprat,[0.10 0.5 0.90]);
LeRBFNet = errorbar(ndata,RBFNetResults(2,:),RBFNetResults(2,:)-RBFNetResults(1,:), RBFNetResults(3,:)-RBFNetResults(2,:),'d--','linewidth',1,'markersize',4);
set(LeRBFNet,'markerfacecolor',get(LeRBFNet,'color'));

legend([LeGP LeAddNet LeFITC LeRBFNet],'Full GP','Additive NN','FITC','RBF NN');
xlabel('Number of observations');
ylabel('Acceptance probability');
xlim([0,3000]);


% Mean ESS per second
figure;
GPResults = quantile(GPMeanESSperSec,[0.10 0.5 0.90]);
LeGP = errorbar(ndata,GPResults(2,:),GPResults(2,:)-GPResults(1,:), GPResults(3,:)-GPResults(2,:),'v--','linewidth',1,'markersize',4);
set(LeGP,'markerfacecolor',get(LeGP,'color'));
hold on;
AddNetResults = quantile(AddNetMeanESSperSec,[0.10 0.5 0.90]);
LeAddNet = errorbar(ndata,AddNetResults(2,:),AddNetResults(2,:)-AddNetResults(1,:), AddNetResults(3,:)-AddNetResults(2,:),'s--','linewidth',1,'markersize',4.5);
set(LeAddNet,'markerfacecolor',get(LeAddNet,'color'));
FITCResults = quantile(FITCMeanESSperSec,[0.10 0.5 0.90]);
LeFITC = errorbar(ndata,FITCResults(2,:),FITCResults(2,:)-FITCResults(1,:), FITCResults(3,:)-FITCResults(2,:),'o--','linewidth',1,'markersize',4);
set(LeFITC,'markerfacecolor',get(LeFITC,'color'));
RBFNetResults = quantile(RBFNetMeanESSperSec,[0.10 0.5 0.90]);
LeRBFNet = errorbar(ndata,RBFNetResults(2,:),RBFNetResults(2,:)-RBFNetResults(1,:), RBFNetResults(3,:)-RBFNetResults(2,:),'d--','linewidth',1,'markersize',4);
set(LeRBFNet,'markerfacecolor',get(LeRBFNet,'color'));

legend([LeGP LeAddNet LeFITC LeRBFNet],'Full GP','Additive NN','FITC','RBF NN');
xlabel('Number of observations');
ylabel('Mean ESS per second');
xlim([0,3000]);

