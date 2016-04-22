function [obj,K,pinK,wts,The] = train(obj,Xtr,Ytr,epsilon)
%if nargin < 4, flag=0;else flag=1; end
if nargin < 4, epsilon = 1e-06; end
%obj.biases{1} = -diag(Xtr*obj.wts{1}');
K = exp(-bsxfun(@rdivide,sq_dist(Xtr',obj.center{1}'),obj.width{1}')/2);
% while det(S) ==0
%     obj = initial(obj,obj.sizes);
%     obj.biases{1} = -diag(Xtr*obj.wts{1}');
%     n = bsxfun(@plus,Xtr*obj.wts{1}',transpose(obj.biases{1}));
%     S = sigmoid(n);
% end
% if flag
%     K = [epsilon*eye(size(K,2));K];
%     Ytr = [zeros(size(K,2),1);Ytr];
% end
% pinK = pinv(K);
% The = pinK*pinK';
%v = pinS*Ytr;
%obj.wts = (pinK*Ytr)';
obj.wts = ((K'*K+epsilon*eye(size(K,2)))\(K'*Ytr))';
%obj.biases{2} = v(1);
wts = obj.wts';