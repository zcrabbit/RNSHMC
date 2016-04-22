function U = BLR_U(y, X, CurrentBeta,alpha,varargin)
if size(CurrentBeta,2)~=1
    CurrentBeta = CurrentBeta';
end
flag = 0;
if nargin > 4
    flag = varargin{1};
    % get subset parameter
    % subset = varargin{2};
end

if flag == 0
    U = sum(log(1+exp(X*CurrentBeta))) - y' * (X*CurrentBeta) + (1/2/alpha)*sum(CurrentBeta.^2);
else
    U = X' *(exp(X*CurrentBeta)./(1+exp(X*CurrentBeta)) - y) + CurrentBeta/alpha;
    % grid parameter space
%     Xsub = X(subset,:);
%     ysub = y(subset);
%     U = Xsub' * (exp(Xsub*CurrentBeta) ./ ( 1+exp(Xsub*CurrentBeta) ) - ysub);
    %U = size(X,1)/length(subset) * U;
    

        
end