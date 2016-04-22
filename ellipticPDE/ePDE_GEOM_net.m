function [dU, G, GDeriv] = ePDE_GEOM_net(net, ProposedTheta,sigmatheta,D)
dU = transpose(gradient(net,ProposedTheta'));
G = (dU-ProposedTheta/sigmatheta^2)*transpose(dU-ProposedTheta/sigmatheta^2)+eye(D)/sigmatheta^2;
d2U = Hessian(net,ProposedTheta)-eye(D)/sigmatheta^2;

for d=1:D
    GDeriv{d} = d2U(:,d)*(dU-ProposedTheta/sigmatheta^2)' + (dU-ProposedTheta/sigmatheta^2)*d2U(d,:);
end