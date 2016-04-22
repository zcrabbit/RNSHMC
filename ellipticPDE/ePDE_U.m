function [output] = ePDE_U(y,theta,PDE,sigmay,sigmatheta, der)
if nargin<6
    der = 0;
end

N=length(y); D=length(theta);
meshsz = sqrt(size(PDE.node,1));
obs_idx=linspace(1,meshsz,sqrt(N));
[obs_mesh1,obs_mesh2]=meshgrid(obs_idx,obs_idx);
idx=sub2ind([meshsz,meshsz],obs_mesh2(:),obs_mesh1(:));
[sol, dsol] = ePDEsolver1(theta,PDE,1);

s = sol(idx); ds = dsol(idx,:);

if der == 0
    loglik = -sum((y-s).^2)/(2*sigmay^2);
    logpri = -theta'*theta/(2*sigmatheta^2);
    output = -(loglik+logpri);
end

if der == 1
    dloglik = sum(repmat(y-s,[1,D]).*ds)'/sigmay^2;
    dlogpri = -theta./sigmatheta^2;
    output = -(dloglik + dlogpri);
end

