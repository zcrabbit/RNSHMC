function [X, Y] = shuffleData(X,Y)
%n = size(X,2);
n = size(X,1);
pi = randperm(n);
% X = X(:,pi);
% Y = Y(:,pi);
X = X(pi,:);
Y = Y(pi,:);