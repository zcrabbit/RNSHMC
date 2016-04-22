function obj = RBFnetworkApprox(sizes,widthscale,center,Sigma)
if nargin < 4, Sigma = eye(sizes(1)); end
% networkApprox(sizes) : construct a neural network based approximation
if nargin < 3, center = zeros(sizes(1),1); end
obj.numoflayers = length(sizes); obj.sizes = sizes;
obj.width = cell(obj.numoflayers-2,1);
obj.center = cell(obj.numoflayers-2,1);
obj.wts = zeros(1,sizes(obj.numoflayers-1));
for numlayer = 1:obj.numoflayers-2
    obj.width{numlayer} = widthscale./chi2rnd(1,sizes(numlayer+1),1);
    %obj.width{numlayer} = widthscale * ones(sizes(numlayer+1),1);
    %obj.center{numlayer} = randn(sizes(numlayer+1),sizes(numlayer));
    obj.center{numlayer} = mvnrnd(zeros(1,sizes(1)),Sigma,sizes(numlayer+1));
end
obj.center{1} = bsxfun(@plus,obj.center{1}*chol(Sigma),center');
%obj.biases{1} = obj.biases{1} - obj.wts{1} * center;
%obj.width{numlayer} = zeros(sizes(numlayer+1),1);
obj=class(obj,'RBFnetworkApprox');
end