function [REM, Times] = ePDE_RMHMC(y,PDE,sigmay,sigmatheta,D,meanTrue)

NumOfIterations = 2000;
BurnIn = floor(0.5*NumOfIterations);
%ColStart= floor(0.1*NumOfIterations);

Trajectory =  2.4;
NumOfLeapFrogSteps = 3;
NumOfNewtonSteps = 5;
StepSize = Trajectory/NumOfLeapFrogSteps;

% H_HMC = zeros((NumOfIterations-BurnIn)*NumOfLeapFrogSteps,1);
thetaSaved = zeros(NumOfIterations-BurnIn,D);

% Initialize
%Startpoint = randn(D,1);
Startpoint = zeros(D,1);
CurrentTheta = Startpoint;
CurrentU = ePDE_U(y,CurrentTheta,PDE,sigmay,sigmatheta);
Times = [0];

% Random numbers
randn('state',2015);
rand('twister',2015);

for d = 1:D
    GDeriv{d} = zeros(D);
end

%% 
Proposed = 0;
Accepted = 0;
thetaSam = [CurrentTheta'];
REM = [norm(CurrentTheta'-meanTrue)/norm(meanTrue)];
tic;

for IterationNum = 1:NumOfIterations
    if mod(IterationNum,100)==0
        disp([num2str(IterationNum) ' iterations completed.'])
        disp(Accepted/Proposed)
        Proposed = 0;
        Accepted = 0;
        drawnow;
    end
   
    
    Proposed = Proposed + 1;
    ProposedTheta = CurrentTheta;
    
    % pre-leapfrog calculation
    % Calculate G and the partial derivatives dG/dw
    [dU,G,GDeriv] = ePDE_GEOM(y,ProposedTheta,PDE,sigmay,sigmatheta, 3);
    %G = diag(diag(G));
    CholG = chol(G);
    InvG  = inv(G);

    for d = 1:D
        InvGdG{d}      = InvG*GDeriv(:,:,d);
        %InvGdG{d}      = InvG*diag(diag(GDeriv(:,:,d)));
        TraceInvGdG(d) = trace(InvGdG{d});
    end
    % terms other than quadratic one
    dphi = dU + 0.5*TraceInvGdG';
    
    % propose momentum
    CurrentMomentum = (randn(1,D)*CholG)';
    ProposedMomentum = CurrentMomentum;
    
    % Calculate current H value
    CurrentLogDet = sum(log(diag(CholG)));
    CurrentH  = CurrentU + CurrentLogDet + (CurrentMomentum'*InvG*CurrentMomentum)/2;
    
    for StepNum = 1:NumOfLeapFrogSteps
        %%%%%%%%%%%%%%%%%%%
        % Update momentum %
        %%%%%%%%%%%%%%%%%%%
        % Multiple fixed point iteration
        PM = ProposedMomentum;
        for FixedIter = 1:NumOfNewtonSteps
            %MomentumHist(FixedIter,:) = PM;
            
            InvGMomentum = InvG*PM;
            for d = 1:D
                dQuadTerm(d)  = 0.5*(PM'*InvGdG{d}*InvGMomentum);
            end
            
            PM = ProposedMomentum + (StepSize/2)*(-dphi + dQuadTerm');
        end
        ProposedMomentum = PM;
        %%%%%%%%%%%%%%%%%%%%%%%
        % Update w parameters %
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Multiple Fixed Point Iteration %%%
        FixedInvGMomentum  = G\ProposedMomentum;
        
        PT = ProposedTheta;
        for FixedIter = 1:NumOfNewtonSteps
            %wHist(FixedIter,:) = PB;
            InvGMomentum = ePDE_Met(y,PT,PDE,sigmay,sigmatheta)\ProposedMomentum;
            %InvGMomentum = diag(diag(ePDE_Met(y,PT,PDE,sigmay,sigmatheta)))\ProposedMomentum;
            PT = ProposedTheta + (StepSize/2)*(FixedInvGMomentum + InvGMomentum);
        end
        ProposedTheta = PT;
        % Update G based on new parameters and the partial derivatives dG/dw
        [dU,G,GDeriv] = ePDE_GEOM(y,ProposedTheta,PDE,sigmay,sigmatheta, 3);
        %G = diag(diag(G));
        InvG = inv(G);

        for d = 1:D
            InvGdG{d}      = InvG*GDeriv(:,:,d);
            %InvGdG{d}      = InvG*diag(diag(GDeriv(:,:,d)));
            TraceInvGdG(d) = trace(InvGdG{d});
        end
        % terms other than quadratic one
        dphi = dU + 0.5*TraceInvGdG';
        %%%%%%%%%%%%%%%%%%%
        % Update momentum %
        %%%%%%%%%%%%%%%%%%%
        InvGMomentum = InvG*ProposedMomentum;
        for d = 1:D
            dQuadTerm(d) = 0.5*(ProposedMomentum'*InvGdG{d}*InvGMomentum);
        end
        ProposedMomentum = ProposedMomentum + (StepSize/2)*(-dphi + dQuadTerm');   

%         ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
%         ProposedTheta = ProposedTheta + StepSize * (InvMass*ProposedMomentum);
%         ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));

        

    end

    
    ProposedMomentum = - ProposedMomentum;
    % calculate potential
    ProposedU = ePDE_U(y,ProposedTheta,PDE,sigmay,sigmatheta);
    
    % Calculate H value
    ProposedLogDet = sum(log(diag(chol(G))));
    ProposedH = ProposedU + ProposedLogDet + (ProposedMomentum'*InvG*ProposedMomentum)/2; 
    
    % calculate the ratio
    Ratio = -ProposedH + CurrentH;
    
    if isfinite(Ratio) && (Ratio > min([0,log(rand)]))
        CurrentTheta = ProposedTheta;
        CurrentU = ProposedU;
        Accepted = Accepted + 1;
    end

    
    % Save samples if required
    if IterationNum > BurnIn
        thetaSaved(IterationNum-BurnIn,:) = CurrentTheta;
    end
    
    thetaSam = [thetaSam; CurrentTheta'];
    REM = [REM norm(mean(thetaSam,1)-meanTrue)/norm(meanTrue)];
    Times = [Times;toc];
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        disp('Burn-in complete, now drawing samples.')
        %tic;
    end
            
end
%Times = toc;

thetaPosterior = thetaSaved;
acprat = size(unique(thetaPosterior,'rows'),1)/(NumOfIterations-BurnIn);

CurTime = fix(clock);
save(['Results/Results_RMHMC_ePDE_'  '_' num2str(CurTime) '.mat'], 'StepSize', 'NumOfLeapFrogSteps', 'acprat', 'thetaPosterior', 'Times')