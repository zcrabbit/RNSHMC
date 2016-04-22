function [A, area] = assembling(node, elem, c)
N = size(node,1); NT = size(elem,1);
ii = zeros(9*NT,1); jj = zeros(9*NT,1); sA = zeros(9*NT,1);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
index = 0;
for i = 1:3
    for j = 1:3
        ii(index+1:index+NT) = elem(:,i);
        jj(index+1:index+NT) = elem(:,j);
        sA(index+1:index+NT) = c.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        index = index + NT;
    end
end
A = sparse(ii,jj,sA,N,N);