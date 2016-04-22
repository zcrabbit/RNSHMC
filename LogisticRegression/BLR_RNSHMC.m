function [betaPosterior,Times,acprat,Trainingdata]= BLR_RNSHMC(s,stepsize,nLeap,data)
% load data
load(['Data/' data]);
% set the regularization for random networks
switch data
    case 'simulation'
        lam = 1e-06;
    case 'a9a60d'
        lam = 1e-09;
    case 'bankmarket'
        lam = 1e-07;
end


NumOfIterations = 10000;
ColStart = 1000;
BurnIn = 5000;

% Trajectory = 0.72;
% NumOfLeapFrogSteps = 20;
% parameters for 'lrdata'
% Trajectory = 1.08;
% NumOfLeapFrogSteps = 24;

% Trajectory = 1.08;
% NumOfLeapFrogSteps = 12;


% parameter for bank marketing data set
% Trajectory = 0.36;
% NumOfLeapFrogSteps = 30;

%parameter for adult data set
%Trajectory = 0.52;
NumOfLeapFrogSteps = nLeap;

alpha = 100;
%StepSize = Trajectory/NumOfLeapFrogSteps;
StepSize = stepsize;
D = size(X,2);
% epsilon for online learning initialization 
% epsilon = 1e-4;


Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));

betaSaved = zeros(NumOfIterations-BurnIn,D);
Xtr = zeros(BurnIn-ColStart,D);
Ytr = zeros(BurnIn-ColStart,1);

Startpoint = zeros(D,1);
% Startpoint_adap = load('Startpoint_lrd28');
% Startpoint_adap = Startpoint_adap.Startpoint;
%Startpoint = [-2;zeros(D-1,1)];
%Startpoint = meanTrue' + 0.5*rand(D,1);
%load Startpoint_lrd28;
CurrentBeta = Startpoint;
CurrentU = BLR_U(y,X,CurrentBeta,alpha);
% Times = [0];

% Random numbers
% randn('state',2015);
% rand('twister',2015);



%%
Proposed = 0;
Accepted = 0;
% ARate = [];
%figure(2); plot(1:length(REM),REM,'b-');
% betaSam = [CurrentBeta'];
% REM = [norm(CurrentBeta'-meanTrue)/norm(meanTrue)];
% tic

for IterationNum = 1:NumOfIterations
        
    if mod(IterationNum,100) == 0
        disp([num2str(IterationNum) ' iterations completed.'])
        disp(Accepted/Proposed)
%         if IterationNum > AdapStart
%             ARate = [ARate Accepted/Proposed];
%             %figure(1); plot((ColStart+100):100:IterationNum,ARate,'*-');
%         end

        Proposed = 0;
        Accepted = 0;
        drawnow
            
    end
    
    ProposedBeta = CurrentBeta;

    % Sample momentum
    ProposedMomentum = (randn(1,D)*Mass)';
    CurrentMomentum = ProposedMomentum;
    
    Proposed = Proposed + 1;

    % Use random Leapfrog steps
    RandomLeapFrogSteps = randsample(NumOfLeapFrogSteps,1);
    % Perform leapfrog steps
    for StepNum = 1:RandomLeapFrogSteps
        if IterationNum <= BurnIn
            ProposedMomentum = ProposedMomentum - StepSize/2.* BLR_U(y,X,ProposedBeta,alpha,1);
            ProposedBeta = ProposedBeta + StepSize.*((InvMass)*ProposedMomentum);
            ProposedMomentum = ProposedMomentum - StepSize/2.*BLR_U(y,X,ProposedBeta,alpha,1);
        else
            ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedBeta'));
            ProposedBeta = ProposedBeta + StepSize * (InvMass*ProposedMomentum);
            ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedBeta'));
        end
    end
    
    ProposedMomentum = -ProposedMomentum;
        
    % Calculate Potential
%     if IterationNum > BurnIn
%         ProposedU = predict(net,ProposedBeta');
%     else
        ProposedU = BLR_U(y,X,ProposedBeta,alpha);
%     end
        
    % Calculate H value
    CurrentH  = CurrentU + (CurrentMomentum'*InvMass*CurrentMomentum)/2;
    ProposedH = ProposedU + (ProposedMomentum'*InvMass*ProposedMomentum)/2;
       
    % Accept according to ratio
    Ratio = -ProposedH + CurrentH;
               
    if (isfinite(Ratio) && (Ratio > min([0,log(rand)])))
        CurrentBeta = ProposedBeta;
        CurrentU = ProposedU;
        Accepted = Accepted + 1;   
    end
    %figure(2); plot(1:length(REM),REM,'b-');

    % Start collection if required
%     if IterationNum > ColStart && IterationNum <= AdapStart
    if IterationNum > ColStart && IterationNum <= BurnIn
        Xtr(IterationNum-ColStart,:) = CurrentBeta;
        Ytr(IterationNum-ColStart) = CurrentU;
    end
    
    % Start training surrogate function
%     if IterationNum == AdapStart
%         %save('Xtrain','Xtr');
%         Xtradp =Xtr(1:(AdapStart-ColStart),:);
%         Ytradp =Ytr(1:(AdapStart-ColStart));
%         base = min(Ytradp); Ytradp = Ytradp - base;
%         net = networkApprox([D,s,1]);
%         %[net,S,pinS] = train(net,Xtradp,Ytradp,epsilon);
%         [net,~,~,v,The] = train(net,Xtradp,Ytradp,epsilon);
% %         CurrentBeta = Startpoint_adap;
% %         CurrentU = BLR_U(y,X,CurrentBeta,alpha);
% %         betaSam = [CurrentBeta'];
% %         REM = [norm(CurrentBeta'-meanTrue)/norm(meanTrue)];
%     end
    
    %if IterationNum > AdapStart && IterationNum <= BurnIn
%     if IterationNum > AdapStart
% %         Xtradp = [Xtradp;CurrentBeta'];
% %         Ytradp = [Ytradp;CurrentU-base];
%         %[S,pinS] = train_online(net,S,pinS,CurrentBeta);
%         [v,The] = train_online_wts(net,v,The,CurrentBeta,CurrentU-base);
%         if rand < (1/sqrt(IterationNum-AdapStart))
%             %net = setwts(net,pinS*[zeros(2000,1);Ytradp]);
%             net = setwts(net,v);
%         end
% %         betaSam = [betaSam; CurrentBeta'];
% %         REM = [REM norm(mean(betaSam,1)-meanTrue)/norm(meanTrue)];
%     end
        
    % Save samples if required
    %betaSaved(IterationNum,:) = CurrentBeta;
    if IterationNum > BurnIn
        betaSaved(IterationNum-BurnIn,:) = CurrentBeta;
    end
    
%     betaSam = [betaSam; CurrentBeta'];
%     REM = [REM norm(mean(betaSam,1)-meanTrue)/norm(meanTrue)];
%     Times = [Times;toc];
    
    % Start timer after burn-in
    if IterationNum == BurnIn
       disp('Burn-in complete, now drawing samples.')
       %CurrentU = predict(net,CurrentBeta');
       base = min(Ytr); Ytr = Ytr - base;
       net = networkApprox([D,s,1]);
%         net = train(net,Xtr,Ytr,0.1);
       net = train(net,Xtr,Ytr,lam);
       tic;
    end
    
end
Times = toc;
Trainingdata.Xtrain = Xtr;
Trainingdata.Ytrain = Ytr;

betaPosterior = betaSaved;
acprat = size(unique(betaPosterior,'rows'),1)/(NumOfIterations-BurnIn);

CurTime = fix(clock);
save(['Results/Results_RNSHMC_BLR' data '_' num2str(nLeap) 'nLeap_' num2str(CurTime) '.mat'], 'StepSize', 'NumOfLeapFrogSteps', 'acprat', 'betaPosterior', 'Times')