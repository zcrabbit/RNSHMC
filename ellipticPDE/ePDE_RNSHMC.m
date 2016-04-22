function [thetaPosterior,Times,acprat,Trainingdata] = ePDE_RNSHMC(s,stepsize,nLeap)
load Data/ePDE
%ColStart= floor(0.001*NumOfIterations);
NumOfIterations = 10000;
BurnIn = 5000;
ColStart = 1000;
% AdapStart = ColStart + ntrain;
% BurnIn = burnin;
% NumOfIterations = BurnIn + 200;

%Trajectory =  2.4;
%NumOfLeapFrogSteps = 10;
NumOfLeapFrogSteps = nLeap;
%StepSize = Trajectory/NumOfLeapFrogSteps;
StepSize = stepsize;
%epsilon = 1e-6;

Mass  = diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));

% H_HMC = zeros((NumOfIterations-BurnIn)*NumOfLeapFrogSteps,1);
thetaSaved = zeros(NumOfIterations-BurnIn,D);
Xtr = zeros(BurnIn-ColStart,D);
Ytr = zeros(BurnIn-ColStart,1);

% Initialize
Startpoint = zeros(D,1);
CurrentTheta = Startpoint;
CurrentU = ePDE_U(y,CurrentTheta,PDE,sigmay,sigmatheta);
%Times = [0];

% % Random numbers
% randn('state',2015);
% rand('twister',2015);

%% 
Proposed = 0;
Accepted = 0;
%ARate = [];
% thetaSam = [CurrentTheta'];
% REM = [norm(CurrentTheta'-meanTrue)/norm(meanTrue)];
% tic;


for IterationNum = 1:NumOfIterations
    if mod(IterationNum,100)==0
        disp([num2str(IterationNum) ' iterations completed.'])
        disp(Accepted/Proposed)
%         if IterationNum > AdapStart
%             ARate = [ARate Accepted/Proposed];
%             %figure(1); plot((ColStart+100):100:IterationNum,ARate,'*-');
%         end
        Proposed = 0;
        Accepted = 0;
        drawnow;
    end
    
    ProposedMomentum = (randn(1,D)*chol(Mass))';
    CurrentMomentum = ProposedMomentum;
    
    Proposed = Proposed + 1;
    ProposedTheta = CurrentTheta;
    
    % Use random Leapfrog steps
    RandomLeapFrogSteps = randsample(NumOfLeapFrogSteps,1);
    % Perform leapfrog steps
    for StepNum = 1:RandomLeapFrogSteps
        if IterationNum <= BurnIn
            ProposedMomentum = ProposedMomentum - StepSize/2.*ePDE_U(y,ProposedTheta,PDE,sigmay,sigmatheta,1);
            ProposedTheta = ProposedTheta + StepSize.*((InvMass)*ProposedMomentum);
            ProposedMomentum = ProposedMomentum - StepSize/2.*ePDE_U(y,ProposedTheta,PDE,sigmay,sigmatheta,1);
        else
            
            ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
            ProposedTheta = ProposedTheta + StepSize * (InvMass*ProposedMomentum);
            ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
        end
           
    end

    
    ProposedMomentum = - ProposedMomentum;
    % calculate potential
    ProposedU = ePDE_U(y,ProposedTheta,PDE,sigmay,sigmatheta);
    
    % calculate Hamiltonian function
    CurrentH = CurrentU + .5*CurrentMomentum'*(InvMass)*CurrentMomentum;
    ProposedH = ProposedU + .5*ProposedMomentum'*(InvMass)*ProposedMomentum;
    
    % calculate the ratio
    Ratio = -ProposedH + CurrentH;
    
    if isfinite(Ratio) && (Ratio > min([0,log(rand)]))
        CurrentTheta = ProposedTheta;
        CurrentU = ProposedU;
        Accepted = Accepted + 1;
    end
    
        % Start collection if required
    if IterationNum > ColStart && IterationNum <= BurnIn
        Xtr(IterationNum-ColStart,:) = CurrentTheta;
        Ytr(IterationNum-ColStart) = CurrentU;
    end
    
%     % Start training surrogate function
%     if IterationNum == AdapStart
%         %save('Xtrain_ePDE','Xtr');
%         Xtradp =Xtr(1:(AdapStart-ColStart),:);
%         Ytradp =Ytr(1:(AdapStart-ColStart));
%         base = min(Ytradp); Ytradp = Ytradp - base;
%         net = networkApprox([D,s,1]);
%         %[net,S,pinS] = train(net,Xtradp,Ytradp,epsilon);
%         [net,~,~,v,The] = train(net,Xtradp,Ytradp,epsilon);
%     end
    
    %if IterationNum > AdapStart && IterationNum <= BurnIn && acceptedornot == 1
%     if IterationNum > AdapStart
% %         Xtradp = [Xtradp;CurrentTheta'];
% %         Ytradp = [Ytradp;CurrentU-base];
%         %[S,pinS] = train_online(net,S,pinS,CurrentBeta);
%         [v,The] = train_online_wts(net,v,The,CurrentTheta,CurrentU-base);
%         if rand < (1/sqrt(IterationNum-AdapStart))
%             %net = setwts(net,pinS*[zeros(2000,1);Ytradp]);
%             %net = train(net,Xtradp,Ytradp,epsilon);
%             net = setwts(net,v);
%         end
%     end

    
    % Save samples if required
    if IterationNum > BurnIn
        thetaSaved(IterationNum-BurnIn,:) = CurrentTheta;
    end
    
%     thetaSam = [thetaSam; CurrentTheta'];
%     REM = [REM norm(mean(thetaSam,1)-meanTrue)/norm(meanTrue)];
%     Times = [Times;toc];
    
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        disp('Burn-in complete, now drawing samples.')
        base = min(Ytr); Ytr = Ytr - base;
        net = networkApprox([D,s,1]);
        %net = train(net,Xtr,Ytr,0.1);
        net = train(net,Xtr,Ytr,1e-07);
        tic;
    end
            
end
Times = toc;
Trainingdata.Xtrain = Xtr;
Trainingdata.Ytrain = Ytr;

thetaPosterior = thetaSaved;
acprat = size(unique(thetaPosterior,'rows'),1)/(NumOfIterations-BurnIn);

CurTime = fix(clock);
save(['Results/Results_RNSHMC_ePDE_'  num2str(nLeap) 'nLeap_' num2str(CurTime) '.mat'], 'StepSize', 'NumOfLeapFrogSteps', 'acprat', 'thetaPosterior', 'Times')