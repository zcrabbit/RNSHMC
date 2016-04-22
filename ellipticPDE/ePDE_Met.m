function [G]=ePDE_Met(y,theta,PDE,sigmay,sigmatheta)

N=length(y); D=length(theta);
meshsz = sqrt(size(PDE.node,1));
obs_idx=linspace(1,meshsz,sqrt(N));
[obs_mesh1,obs_mesh2]=meshgrid(obs_idx,obs_idx);
idx=sub2ind([meshsz,meshsz],obs_mesh2(:),obs_mesh1(:));
[sol, dsol] = ePDEsolver1(theta,PDE,1);

s = sol(idx); ds = dsol(idx,:);


G = (ds'*ds)/sigmay^2 + eye(D)/sigmatheta^2;
