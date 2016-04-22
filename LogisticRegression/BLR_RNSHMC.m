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

NumOfLeapFrogSteps = nLeap;

alpha = 100;
StepSize = stepsize;
D = size(X,2);

Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));

betaSaved = zeros(NumOfIterations-BurnIn,D);
Xtr = zeros(BurnIn-ColStart,D);
Ytr = zeros(BurnIn-ColStart,1);

Startpoint = zeros(D,1);
CurrentBeta = Startpoint;
CurrentU = BLR_U(y,X,CurrentBeta,alpha);

Proposed = 0;
Accepted = 0;

for IterationNum = 1:NumOfIterations
        
    if mod(IterationNum,100) == 0
        disp([num2str(IterationNum) ' iterations completed.'])
        disp(Accepted/Proposed)

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

    % Start collection if required
    if IterationNum > ColStart && IterationNum <= BurnIn
        Xtr(IterationNum-ColStart,:) = CurrentBeta;
        Ytr(IterationNum-ColStart) = CurrentU;
    end
        
    % Save samples if required
    if IterationNum > BurnIn
        betaSaved(IterationNum-BurnIn,:) = CurrentBeta;
    end
    
    % Start timer after burn-in
    if IterationNum == BurnIn
       disp('Burn-in complete, now drawing samples.')
       base = min(Ytr); Ytr = Ytr - base;
       net = networkApprox([D,s,1]);
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
