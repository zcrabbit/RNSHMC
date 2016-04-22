function [thetasaved,Times,Data_training]=Banana_HMC(y)
NumOfIterations = 1000;
BurnIn = floor(0.5*NumOfIterations);
ColStart = floor(0.10*NumOfIterations);

% Random numbers
rng(2016)

% HMC parameters
Trajectory = 1.45;
NumOfLeapFrogSteps = 15;

StepSize = Trajectory/NumOfLeapFrogSteps;
D = 2;

Xtr = zeros(BurnIn-ColStart,D);
Ytr = zeros(BurnIn-ColStart,1);

Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));
thetasaved = zeros(NumOfIterations-BurnIn,D);
CurrentTheta = [-1;1.5];
%EarlyTra = CurrentTheta';
CurrentU = U(y,CurrentTheta);
Accepted = 0;
Proposed = 0;

for iter = 1:NumOfIterations
    % Sample Momentum
    ProposedTheta = CurrentTheta;
%     if iter <= 10
%         EarlyTra = [EarlyTra;CurrentTheta'];
%     end
    CurrentMomentum = (randn(1,D)*Mass)';
    ProposedMomentum = CurrentMomentum;
    
    
    % Perform Leapfrog Steps
    for StepNum = 1:NumOfLeapFrogSteps
%         if iter <= BurnIn
        ProposedMomentum = ProposedMomentum - StepSize/2 * U(y,ProposedTheta,1); 
        ProposedTheta = ProposedTheta + StepSize * (InvMass*ProposedMomentum);
        ProposedMomentum = ProposedMomentum - StepSize/2 * U(y,ProposedTheta,1);
%         else
%             ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
%             ProposedTheta = ProposedTheta + StepSize * (InvMass*ProposedMomentum);
%             ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
%         end
%         
%         if iter <=20
%             EarlyTra = [EarlyTra;ProposedTheta'];
%         end
    end
    
    ProposedMomentum = - ProposedMomentum;
    ProposedU = U(y,ProposedTheta);
    
    % Calculate Hamiltonians
    CurrentH = CurrentU + .5*CurrentMomentum'*InvMass*CurrentMomentum;
    ProposedH = ProposedU + .5*ProposedMomentum'*InvMass*ProposedMomentum;
    
    % Metropolis-Hasting step
    ratio = CurrentH - ProposedH;
    if isfinite(ratio) && ratio > min([0,log(rand)])
        CurrentTheta = ProposedTheta;
        CurrentU = ProposedU;
        Accepted = Accepted + 1;
%     else if iter <= 20
%         EarlyTra((end-NumOfLeapFrogSteps+1):end,:) = [];
%         end
    end
    
    % Start collection if required
    if iter > ColStart && iter <= BurnIn
        Xtr(iter-ColStart,:) = CurrentTheta;
        Ytr(iter-ColStart) = CurrentU;
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
%         base = min(Ytr); Ytr = Ytr - base;
%         net = networkApprox([D,40,1]);
% %         net = train(net,Xtr,Ytr,0.1);
%         net = train(net,Xtr,Ytr);
        tic;
    end
    
end     
Times = toc;
Data_training.Xtr = Xtr;
Data_training.Ytr = Ytr;

% markerind = 1:NumOfLeapFrogSteps:size(EarlyTra,1);
% xmarkers = EarlyTra(markerind,1);
% ymarkers = EarlyTra(markerind,2);
% 
% plot(EarlyTra(:,1),EarlyTra(:,2),'-r',xmarkers,ymarkers,'ob',...
%                         'MarkerSize',2,'MarkerEdgeColor','b',...
%                         'MarkerFaceColor','b');
end


function z = U(y,theta,der)
if nargin <=2, der =0;end
n = size(y,1);
if der ==0
    z = 0.5*sum((y-theta(1)-theta(2)^2).^2)/4 + 0.5*sum(theta.^2);  
else
    z = (n*(theta(1)+theta(2)^2)-sum(y))/4 *[1;2*theta(2)] + theta;
end
end
