function [thetasaved,Time,acprat,TrainTime] = NarrowGaussian_HMC(Sigma,ndata,method)
NumOfIterations = 2000;
BurnIn = floor(0.2*NumOfIterations);
D = size(Sigma,1);
nu = 1000;

% Random numbers
% randn('state',2015);
% rand('twister',2015);
rng(2015);

% training data
Xtrain = mvnrnd(zeros(1,32), Sigma,ndata);
Ytrain = .5*sum((Xtrain/Sigma) .* Xtrain,2)+0.1*randn(size(Xtrain,1),1);

% training surrogate
switch method
    case 'GP'
        hyp.cov = [0,0]; hyp.lik = log(0.1); hyp.mean = [];
        covfunc = @covSEiso; 
        likfunc = @likGauss;
        hyp = minimize(hyp, @gp, -1000, @infExact, [], covfunc, likfunc, Xtrain, Ytrain);
        tic;
        post = infExact(hyp,{@meanZero},{covfunc},likfunc,Xtrain,Ytrain);
        alpha = post.alpha;
        TrainTime = toc;
    case 'NN'
        tic;
        base = min(Ytrain); Ytrain = Ytrain - base;
        net = networkApprox([D,nu,1]);
        net = train(net,Xtrain,Ytrain);
        TrainTime = toc;
    
    case 'FITC'
        hyp2.cov = [0,0]; hyp2.lik = log(0.1); hyp2.mean = [];
        covfunc = @covSEiso; 
        likfunc = @likGauss;
        Xu = mvnrnd(zeros(1,D), Sigma,1000);
        covfuncF = {@covFITC, {covfunc}, Xu};
        hyp2 = minimize(hyp2, @gp, -1000, @infFITC, [], covfuncF, likfunc, Xtrain, Ytrain);
        tic;
        post = infFITC(hyp2,{@meanZero},covfuncF,likfunc,Xtrain,Ytrain);
        alpha = post.alpha;
        TrainTime = toc;
        
    case 'RBF'
        tic;
        net = RBFnetworkApprox([D,1000,1],1,zeros(32,1),Sigma);
        net = train(net,Xtrain,Ytrain);
        TrainTime = toc;
end
        
       
% HMC parameters
NumOfLeapFrogSteps = 30;

StepSize = 0.108;


Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));
thetasaved = zeros(NumOfIterations-BurnIn,D);
CurrentTheta = zeros(1,D);
CurrentU = .5*CurrentTheta*(Sigma\CurrentTheta');
Accepted = 0;
Proposed = 0;

rng('shuffle');

for iter = 1:NumOfIterations
    % Sample Momentum
    ProposedTheta = CurrentTheta;
    CurrentMomentum = (randn(1,D)*Mass);
    ProposedMomentum = CurrentMomentum;
    
    
    % Perform Leapfrog Steps
    for StepNum = 1:NumOfLeapFrogSteps
        switch method
            case 'Exact'
                ProposedMomentum = ProposedMomentum - StepSize/2 * (ProposedTheta/Sigma);
                ProposedTheta = ProposedTheta + StepSize * (ProposedMomentum*InvMass);
                ProposedMomentum = ProposedMomentum - StepSize/2 * (ProposedTheta/Sigma);
                
            case 'GP'
                ProposedMomentum = ProposedMomentum - StepSize/2 * gp_grad(hyp,covfunc,alpha,Xtrain,ProposedTheta);
                ProposedTheta = ProposedTheta + StepSize * (ProposedMomentum*InvMass);
                ProposedMomentum = ProposedMomentum - StepSize/2 * gp_grad(hyp,covfunc,alpha,Xtrain,ProposedTheta);
                
            case 'NN'
                ProposedMomentum = ProposedMomentum - StepSize/2 * (gradient(net,ProposedTheta));
                ProposedTheta = ProposedTheta + StepSize * (ProposedMomentum*InvMass);
                ProposedMomentum = ProposedMomentum - StepSize/2 * (gradient(net,ProposedTheta));
            
            case 'FITC'
                ProposedMomentum = ProposedMomentum - StepSize/2 * gp_grad(hyp2,covfuncF,alpha,Xtrain,ProposedTheta);
                ProposedTheta = ProposedTheta + StepSize * (ProposedMomentum*InvMass);
                ProposedMomentum = ProposedMomentum - StepSize/2 * gp_grad(hyp2,covfuncF,alpha,Xtrain,ProposedTheta);
                
            case 'RBF'
                ProposedMomentum = ProposedMomentum - StepSize/2 * (gradient(net,ProposedTheta));
                ProposedTheta = ProposedTheta + StepSize * (ProposedMomentum*InvMass);
                ProposedMomentum = ProposedMomentum - StepSize/2 * (gradient(net,ProposedTheta));
                
        end
    end
    
    ProposedMomentum = - ProposedMomentum;
    ProposedU = .5*ProposedTheta*(Sigma\ProposedTheta');
    
    % Calculate Hamiltonians
    CurrentH = CurrentU + .5*CurrentMomentum*InvMass*CurrentMomentum';
    ProposedH = ProposedU + .5*ProposedMomentum*InvMass*ProposedMomentum';
    
    % Metropolis-Hasting step
    ratio = CurrentH - ProposedH;
    if isfinite(ratio) && ratio > min([0,log(rand)])
        CurrentTheta = ProposedTheta;
        CurrentU = ProposedU;
        Accepted = Accepted + 1;
    end
    
    Proposed = Proposed + 1;
        
    
    if mod(iter,100)==0
        disp([num2str(iter) ' iterations completed.'])
        disp(Accepted/Proposed)
        
        Accepted = 0;
        Proposed = 0;
    end
    
    if iter > BurnIn
        thetasaved(iter-BurnIn,:) = CurrentTheta;
    end
    
    
    % Start timer after burn-in
    if iter == BurnIn
        disp('Burn-in complete, now drawing samples.')
        tic;
    end
    
end     
Time = toc;
acprat = size(unique(thetasaved,'rows'),1)/(NumOfIterations-BurnIn);
end
