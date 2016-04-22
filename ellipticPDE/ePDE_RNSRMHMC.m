function [thetaPosterior,Times,acprat] = ePDE_RNSRMHMC(stepsize,nLeap)
%%
% load data
load Data/ePDE
NumOfIterations = 10000;
BurnIn = floor(0.5*NumOfIterations);
ColStart= floor(0.1*NumOfIterations);

NumOfLeapFrogSteps = nLeap;
NumOfNewtonSteps = 5;
StepSize = stepsize;

thetaSaved = zeros(NumOfIterations-BurnIn,D);
Xtr = zeros(BurnIn-ColStart,D);
Ytr = zeros(BurnIn-ColStart,1);

% Initialize
Startpoint = randn(D,1);
CurrentTheta = Startpoint;
CurrentU = ePDE_U(y,CurrentTheta,PDE,sigmay,sigmatheta);

for d = 1:D
    GDeriv{d} = zeros(D);
end

%% 
Proposed = 0;
Accepted = 0;

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
    if IterationNum <= BurnIn
        [dU,G,GDeriv] = ePDE_GEOM(y,ProposedTheta,PDE,sigmay,sigmatheta, 3);
    else
        dU = transpose(gradient(net,ProposedTheta'));
        G = diag(diag(Hessian(net,ProposedTheta)));
        GDeriv = Hessian(net,ProposedTheta,1);
        
    end
    CholG = chol(G);
    InvG  = inv(G);

    for d = 1:D
        if IterationNum <= BurnIn
            InvGdG{d} = InvG*GDeriv(:,:,d);
        else
            InvGdG{d} = InvG*diag(diag(GDeriv{d}));
        end
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
            if IterationNum <= BurnIn
                InvGMomentum = ePDE_Met(y,PT,PDE,sigmay,sigmatheta)\ProposedMomentum;
            else
                InvGMomentum = diag(diag(Hessian(net,PT)))\ProposedMomentum;
            end
            PT = ProposedTheta + (StepSize/2)*(FixedInvGMomentum + InvGMomentum);
        end
        ProposedTheta = PT;
        % Update G based on new parameters and the partial derivatives dG/dw
        if IterationNum <= BurnIn
            [dU,G,GDeriv] = ePDE_GEOM(y,ProposedTheta,PDE,sigmay,sigmatheta, 3);
        else
            G = diag(diag(Hessian(net,ProposedTheta)));
            GDeriv = Hessian(net,ProposedTheta,1);
            dU = transpose(gradient(net,ProposedTheta'));
        end
            
        InvG = inv(G);

        for d = 1:D
            if IterationNum <= BurnIn
                InvGdG{d} = InvG*GDeriv(:,:,d);
            else
                InvGdG{d} = InvG*diag(diag(GDeriv{d}));
            end
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
    
    % Start collection if required
    if IterationNum > ColStart && IterationNum <= BurnIn
        Xtr(IterationNum-ColStart,:) = CurrentTheta;
        Ytr(IterationNum-ColStart) = CurrentU;
    end
   
    % Save samples if required
    if IterationNum > BurnIn
        thetaSaved(IterationNum-BurnIn,:) = CurrentTheta;
    end
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        disp('Burn-in complete, now drawing samples.')
        base = min(Ytr); Ytr = Ytr - base;
        net = networkApprox([D,1000,1]);
        net = train(net,Xtr,Ytr,1e-02);
        tic;
    end
            
end

%%
Times = toc;

thetaPosterior = thetaSaved;
acprat = size(unique(thetaPosterior,'rows'),1)/(NumOfIterations-BurnIn);

CurTime = fix(clock);
save(['Results/Results_RNSRMHMC_ePDE_' num2str(nLeap)  'nLeap_' num2str(CurTime) '.mat'], 'StepSize', 'NumOfLeapFrogSteps', 'acprat', 'thetaPosterior', 'Times')