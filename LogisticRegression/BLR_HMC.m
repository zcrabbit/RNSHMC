function [betaPosterior,Times,acprat] = BLR_HMC(stepsize,nLeaps,data)
% load data
load(['Data/' data]);
NumOfIterations = 5000;
BurnIn = floor(0.2*NumOfIterations);


% Trajectory = 0.52;
% NumOfLeapFrogSteps = 40;

% parameters for 'lrdata'
% Trajectory = 1.08;
% NumOfLeapFrogSteps = 24;

% Trajectory = 1.08;
% NumOfLeapFrogSteps = 12;

% Trajectory = 0.33;
% NumOfLeapFrogSteps = 30;

%parameter for adult data set
%Trajectory = 0.52;
NumOfLeapFrogSteps = nLeaps;

alpha = 100;
%StepSize = Trajectory/NumOfLeapFrogSteps;
StepSize = stepsize;
D = size(X,2);


Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));

betaSaved = zeros(NumOfIterations-BurnIn,D);
% Xtr = zeros(BurnIn-ColStart,D);
% Ytr = zeros(BurnIn-ColStart,1);

Startpoint = zeros(D,1);
% Startpoint_adap = load('Startpoint_lrd28');
% Startpoint_adap = Startpoint_adap.Startpoint;
%Startpoint = [-2;zeros(D-1,1)];
%Startpoint = [1;-1;zeros(D-2,1)];
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
% betaSam = [CurrentBeta'];
% REM = [norm(CurrentBeta'-meanTrue)/norm(meanTrue)];
% tic


for IterationNum = 1:NumOfIterations   
    if mod(IterationNum,100) == 0
        disp([num2str(IterationNum) ' iterations completed.'])
        disp(Accepted/Proposed)
%         if IterationNum > BurnIn
%             ARate = [ARate Accepted/Proposed];
%             %figure(1); plot((ColStart+100):100:IterationNum,ARate,'*-');
%             %figure(1); semilogy(REM);
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
        ProposedMomentum = ProposedMomentum - StepSize/2.*BLR_U(y,X,ProposedBeta,alpha,1);
        ProposedBeta = ProposedBeta + StepSize.*((InvMass)*ProposedMomentum);
        ProposedMomentum = ProposedMomentum - StepSize/2.*BLR_U(y,X,ProposedBeta,alpha,1);
    end
    
    ProposedMomentum = -ProposedMomentum;
        
    % Calculate Potential
    ProposedU = BLR_U(y,X,ProposedBeta,alpha);
        
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
    
%     betaSam = [betaSam; CurrentBeta'];
%     REM = [REM norm(mean(betaSam,1)-meanTrue)/norm(meanTrue)];
%     Times = [Times;toc];

    % Start collection if required
%     if IterationNum > ColStart && IterationNum <= BurnIn
%         Xtr(IterationNum-ColStart,:) = CurrentBeta;
%         Ytr(IterationNum-ColStart) = CurrentU;
%     end
        
    % Save samples if required
    %betaSaved(IterationNum,:) = CurrentBeta;
    if IterationNum > BurnIn
        betaSaved(IterationNum-BurnIn,:) = CurrentBeta;
    end
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        disp('Burn-in complete, now drawing samples.')
%         CurrentBeta = Startpoint_adap;
%         CurrentU = BLR_U(y,X,CurrentBeta,alpha);
%         betaSam = [CurrentBeta'];
%         REM = [norm(CurrentBeta'-meanTrue)/norm(meanTrue)];
%         base = min(Ytr); Ytr = Ytr - base;
%         net = networkApprox([50,2000,1]);
%         net = train(net,Xtr,Ytr,0.1);
        tic;
    end
end
Times = toc;

betaPosterior = betaSaved;
acprat = size(unique(betaPosterior,'rows'),1)/(NumOfIterations-BurnIn);

CurTime = fix(clock);
save(['Results/Results_HMC_BLR' data  '_' num2str(NumOfLeapFrogSteps) 'nLeap_' num2str(stepsize) 'stepsize_' num2str(CurTime) '.mat'], 'StepSize', 'NumOfLeapFrogSteps', 'acprat', 'betaPosterior', 'Times')