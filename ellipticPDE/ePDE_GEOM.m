function [dU,G,dG]=ePDE_GEOM(y,theta,PDE,sigmay,sigmatheta, der)
if nargin<6
    der = 0;
end

N=length(y); D=length(theta);
meshsz = sqrt(size(PDE.node,1));
obs_idx=linspace(1,meshsz,sqrt(N));
[obs_mesh1,obs_mesh2]=meshgrid(obs_idx,obs_idx);
idx=sub2ind([meshsz,meshsz],obs_mesh2(:),obs_mesh1(:));
[sol, dsol,d2sol] = ePDEsolver1(theta,PDE,2);

s = sol(idx); ds = dsol(idx,:); d2s = d2sol(idx,:,:);

if der >= 1
    dloglik = sum(repmat(y-s,[1,D]).*ds)'/sigmay^2;
    dlogpri = -theta./sigmatheta^2;
    dU = -(dloglik + dlogpri);
    if der >= 2
        G = (ds'*ds)/sigmay^2 + eye(D)/sigmatheta^2;
        if der == 3
            dG = reshape(ds'*d2s(:,:),repmat(D,1,3));
            dG = dG + permute(dG,[2,1,3]);
            dG = dG./sigmay^2;
        end
    end
    
end