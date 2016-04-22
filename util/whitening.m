function X = whitening(X)
sigma = X'*X/size(X,1);
[U,S,V] = svd(sigma);
X = X * U * diag(1./sqrt(diag(S)+1e-9));