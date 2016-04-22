%% Compare the trajectories of HMC and RNSHMC on a 2D Banana-shaped distribution
%% generate data
rng(2016)
y = 1 + 2*randn(100,1);

%% contour plot of the target density function
TargetDen = @(theta) Banana_U(y,theta);
Dom = [-2 2];
figure(1);
contourplot(TargetDen,Dom);

%% force plot of the target
Force = @(theta) Banana_U(y,theta,1);
hold on;
forceplot(Force,Dom,'r');

%% Run HMC to collect training data and train the random network surrogate
[thetasaved,Times,Data_training] = Banana_HMC(y);
Xtrain = Data_training.Xtr; Ytrain = Data_training.Ytr;
base = min(Ytrain); Ytrain = Ytrain - base;
net = networkApprox([2,40,1]);
net = train(net,Xtrain,Ytrain);

%% compare the force maps
RNForce = @(theta) gradient(net,theta);
hold on;
forceplot(RNForce,Dom,'b');

%% compare the trajectories
% initialize momentums
nTraj = 5;
Moment = randn(2,nTraj);

% HMC
figure(2);
contourplot(TargetDen,Dom);
hold on;
Banana_traj(y,nTraj,Moment);

%RNSHMC
figure(3);
contourplot(TargetDen,Dom);
hold on;
Banana_traj(y,nTraj,Moment,net);

