function [output]=BLR_Met(XX,w,alpha,der)
if(nargin<4)
    der=0;
end

[N,D] = size(XX);
f = XX*w;
p = 1./(1+exp(-f));
v = p.*(ones(N,1)-p);

% Calculate G
if der==0
    G = (XX'.*repmat(v',D,1))*XX + (eye(D)./alpha);
    output = G;
% Calculate the partial derivatives dG/dw
% Faster to compute Z as vector then diag in the next calculation
elseif der==1
    for d = 1:D
        Z = ((1 - 2*p).*XX(:,d));
        %GDeriv{d} = XX'*v*diag(Z)*XX; % Slow because of diag
        % faster
        Z1 = (v.*Z);
        for a =1:D Z2(:,a) = (XX(:,a).*Z1); end
        GDeriv{d} = Z2'*XX;
    end
    output = GDeriv;
else
    disp('wrong choice of der!')
end