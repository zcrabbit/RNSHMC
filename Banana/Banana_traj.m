function Banana_traj(y,T,Moment,net)
% parameters
NumOfIterations = T;
Trajectory = 1.45;
NumOfLeapFrogSteps = 40;

StepSize = Trajectory/NumOfLeapFrogSteps;
D = 2;

Mass  =  diag(ones(D,1)*1);
InvMass = sparse(inv(Mass));
CurrentTheta = [-1;1.5];


for iter = 1:NumOfIterations
    % Sample Momentum
    ProposedTheta = CurrentTheta;
    Traj = CurrentTheta';
    CurrentMomentum = Moment(:,iter);
    ProposedMomentum = CurrentMomentum;
    
    
    % Perform Leapfrog Steps
    for StepNum = 1:NumOfLeapFrogSteps
        if nargin <= 3
            ProposedMomentum = ProposedMomentum - StepSize/2 * U(y,ProposedTheta,1);
            ProposedTheta = ProposedTheta + StepSize * (InvMass*ProposedMomentum);
            ProposedMomentum = ProposedMomentum - StepSize/2 * U(y,ProposedTheta,1);
        else
            ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
            ProposedTheta = ProposedTheta + StepSize * (InvMass*ProposedMomentum);
            ProposedMomentum = ProposedMomentum - StepSize/2 * transpose(gradient(net,ProposedTheta'));
        end
        
        Traj = [Traj;ProposedTheta'];
    end
    
    markerind = 1:NumOfLeapFrogSteps:size(Traj,1);
    xmarkers = Traj(markerind,1);
    ymarkers = Traj(markerind,2);
    
    plot(Traj(:,1),Traj(:,2),'-r',xmarkers(1),ymarkers(1),'pr',...
        'MarkerSize',10,'MarkerEdgeColor','r',...
        'MarkerFaceColor','r'); hold on;
    plot(xmarkers(2),ymarkers(2),'ob','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b');
    Traj = [];
    
end
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