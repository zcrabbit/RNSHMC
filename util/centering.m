function X = centering(X)
center = mean(X,1);
X = bsxfun(@minus,X,center);