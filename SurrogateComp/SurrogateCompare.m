%% Compare surrogate functions on 1D Gaussian distribution.
% set random seed
rng(12345);
n = 5;
%% Case 1: N = 10
%X = randn(10,1);
load Data/dataN10
% Y = X.^2/2 + 0.1*randn(size(X));
% save('dataN10t','X','Y');

figure(1); plot(X,Y,'.','markersize',15);
z = linspace(-4,4,201)';
hold on; 
LeTrue = plot(z,z.^2/2,'linewidth',1);

% Full GP Surrogate
hyp.cov = [0,0]; hyp.lik = log(0.1);hyp.mean = [];
covfunc = @covSEiso;
likfunc = @likGauss;
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X, Y);
m = gp(hyp,@infExact,[],covfunc,likfunc,X,Y,z);
LeGP = plot(z,m,'linewidth',1);

% NN with additive nodes
mAddNet = zeros(size(z,1),20);
for i = 1:100
    netAdd = networkApprox([1,n,1]);
    netAdd = train(netAdd,X,Y);
    mAddNet(:,i) = predict(netAdd,z);
end
LeAddNet = plot(z,mean(mAddNet,2),'linewidth',1);

% FITC surrogate
mF = zeros(size(z,1),20);
for i = 1:100
    hyp2.cov = [0,0];hyp2.lik = log(0.1);hyp.mean = [];
    % Assign inducing points
    Xu = randn(n,1);
    covfuncF = {@covFITC, {covfunc}, Xu};
    hyp2 = minimize(hyp2, @gp, -100, @infFITC, [], covfuncF, likfunc, X, Y);
    mF(:,i) = gp(hyp2,@infFITC,[],covfuncF,likfunc,X,Y,z);
end
LeTICF = plot(z,mean(mF,2),'linewidth',1);

% NN with radial basis function (RBF)
mRBFNet = zeros(size(z,1),20);
for i = 1:100
    netRBF = RBFnetworkApprox([1,n,1],exp(2*hyp2.cov(1)));
    netRBF = train(netRBF,X,Y);
    mRBFNet(:,i) = predict(netRBF,z);
end
LeRBFNet = plot(z,mean(mRBFNet,2),'linewidth',1);

legend([LeTrue LeGP LeAddNet LeTICF LeRBFNet],'Exact','Full GP','Additive NN','TICF','RBF NN');


%% Case 2: N = 20
load Data/dataN20
% Xnew = randn(10,1);
% Ynew = Xnew.^2/2 + 0.1*randn(size(Xnew));
% X = [X;Xnew];
% Y = [Y;Ynew];
% save('dataN20t','X','Y');

figure(2); plot(X,Y,'.','markersize',15);
z = linspace(-4,4,201)';
hold on; 
LeTrue = plot(z,z.^2/2,'linewidth',1);

% Full GP Surrogate
hyp.cov = [0,0]; hyp.lik = log(0.1);hyp.mean = [];
covfunc = @covSEiso;
likfunc = @likGauss;
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X, Y);
m = gp(hyp,@infExact,[],covfunc,likfunc,X,Y,z);
LeGP = plot(z,m,'linewidth',1);

% NN with additive nodes
mAddNet = zeros(size(z,1),20);
for i = 1:100
    netAdd = networkApprox([1,n,1]);
    netAdd = train(netAdd,X,Y);
    mAddNet(:,i) = predict(netAdd,z);
end
LeAddNet = plot(z,mean(mAddNet,2),'linewidth',1);

% FITC surrogate
mF = zeros(size(z,1),20);
for i = 1:100
    hyp2.cov = [0,0];hyp2.lik = log(0.1);hyp.mean = [];
    % Assign inducing points
    Xu = randn(n,1);
    covfuncF = {@covFITC, {covfunc}, Xu};
    hyp2 = minimize(hyp2, @gp, -100, @infFITC, [], covfuncF, likfunc, X, Y);
    mF(:,i) = gp(hyp2,@infFITC,[],covfuncF,likfunc,X,Y,z);
end
LeTICF = plot(z,mean(mF,2),'linewidth',1);

% NN with radial basis function (RBF)
mRBFNet = zeros(size(z,1),20);
for i = 1:100
    netRBF = RBFnetworkApprox([1,n,1],exp(2*hyp2.cov(1)));
    netRBF = train(netRBF,X,Y);
    mRBFNet(:,i) = predict(netRBF,z);
end
LeRBFNet = plot(z,mean(mRBFNet,2),'linewidth',1);

legend([LeTrue LeGP LeAddNet LeTICF LeRBFNet],'Exact','Full GP','Additive NN','TICF','RBF NN');


%% Case 3: N = 40
load Data/dataN40
% Xnew = randn(20,1);
% Ynew = Xnew.^2/2 + 0.1*randn(size(Xnew));
% X = [X;Xnew];
% Y = [Y;Ynew];

figure(3); plot(X,Y,'.','markersize',15);
z = linspace(-4,4,201)';
hold on; 
LeTrue = plot(z,z.^2/2,'linewidth',1);

% Full GP Surrogate
hyp.cov = [0,0]; hyp.lik = log(0.1);hyp.mean = [];
covfunc = @covSEiso;
likfunc = @likGauss;
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X, Y);
m = gp(hyp,@infExact,[],covfunc,likfunc,X,Y,z);
LeGP = plot(z,m,'linewidth',1);

% NN with additive nodes
mAddNet = zeros(size(z,1),20);
for i = 1:100
    netAdd = networkApprox([1,n,1]);
    netAdd = train(netAdd,X,Y);
    mAddNet(:,i) = predict(netAdd,z);
end
LeAddNet = plot(z,mean(mAddNet,2),'linewidth',1);

% FITC surrogate
mF = zeros(size(z,1),20);
for i = 1:100
    hyp2.cov = [0,0];hyp2.lik = log(0.1);hyp.mean = [];
    % Assign inducing points
    Xu = randn(n,1);
    covfuncF = {@covFITC, {covfunc}, Xu};
    hyp2 = minimize(hyp2, @gp, -100, @infFITC, [], covfuncF, likfunc, X, Y);
    mF(:,i) = gp(hyp2,@infFITC,[],covfuncF,likfunc,X,Y,z);
end
LeTICF = plot(z,mean(mF,2),'linewidth',1);

% NN with radial basis function (RBF)
mRBFNet = zeros(size(z,1),20);
for i = 1:100
    netRBF = RBFnetworkApprox([1,n,1],exp(2*hyp2.cov(1)));
    netRBF = train(netRBF,X,Y);
    mRBFNet(:,i) = predict(netRBF,z);
end
LeRBFNet = plot(z,mean(mRBFNet,2),'linewidth',1);

legend([LeTrue LeGP LeAddNet LeTICF LeRBFNet],'Exact','Full GP','Additive NN','TICF','RBF NN');