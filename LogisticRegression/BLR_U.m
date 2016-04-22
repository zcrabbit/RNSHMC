function U = BLR_U(y, X, CurrentBeta,alpha,varargin)
if size(CurrentBeta,2)~=1
    CurrentBeta = CurrentBeta';
end
flag = 0;
if nargin > 4
    flag = varargin{1};
end

if flag == 0
    U = sum(log(1+exp(X*CurrentBeta))) - y' * (X*CurrentBeta) + (1/2/alpha)*sum(CurrentBeta.^2);
else
    U = X' *(exp(X*CurrentBeta)./(1+exp(X*CurrentBeta)) - y) + CurrentBeta/alpha;

end
