function [betaPosterior,Times,acprat] = BLR_HMC(stepsize,nLeaps,data)
% load data
load(['Data/' data]);
NumOfIterations = 5000;
BurnIn = floor(0.2*NumOfIterations);

NumOfLeapFrogSteps = nLeaps;

alpha = 100;
StepSize = stepsize;
D = size(X,2);


Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));

betaSaved = zeros(NumOfIterations-BurnIn,D);

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
        
    % Save samples if required
    if IterationNum > BurnIn
        betaSaved(IterationNum-BurnIn,:) = CurrentBeta;
    end
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        disp('Burn-in complete, now drawing samples.')
        tic;
    end
end
Times = toc;

betaPosterior = betaSaved;
acprat = size(unique(betaPosterior,'rows'),1)/(NumOfIterations-BurnIn);

CurTime = fix(clock);
save(['Results/Results_HMC_BLR' data  '_' num2str(NumOfLeapFrogSteps) 'nLeap_' num2str(stepsize) 'stepsize_' num2str(CurTime) '.mat'], 'StepSize', 'NumOfLeapFrogSteps', 'acprat', 'betaPosterior', 'Times')