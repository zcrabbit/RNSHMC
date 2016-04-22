function obj = networkApprox(sizes,center,Sigma)
if nargin < 3, Sigma = eye(sizes(1)); end
% networkApprox(sizes) : construct a neural network based approximation
if nargin < 2, center = zeros(sizes(1),1); end
obj.numoflayers = length(sizes); obj.sizes = sizes;
obj.biases = cell(obj.numoflayers-1,1);
obj.wts = cell(obj.numoflayers-1,1);
for numlayer = 1:obj.numoflayers-1
    obj.biases{numlayer} = randn(sizes(numlayer+1),1);
    obj.wts{numlayer} = randn(sizes(numlayer+1),sizes(numlayer))/sqrt(sizes(numlayer));
end
obj.wts{1} = obj.wts{1}*chol(Sigma);
obj.biases{1} = obj.biases{1} - obj.wts{1} * center;
obj.biases{numlayer} = zeros(sizes(numlayer+1),1);
obj=class(obj,'networkApprox');
end