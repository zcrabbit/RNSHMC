function U = Banana_U(y,theta,der)
if nargin <=2, der =0;end
n = size(y,1);
m = size(theta,1);
mu = sum([theta(:,1) theta(:,2).^2],2);
if der == 0
    U = 0.5*sum((repmat(y',m,1)-repmat(mu,1,n)).^2,2)/4 + 0.5*sum(theta.^2,2);
else
    U = bsxfun(@times,(n*mu-sum(y))/4, [ones(m,1) 2*theta(:,2)]) + theta;
end