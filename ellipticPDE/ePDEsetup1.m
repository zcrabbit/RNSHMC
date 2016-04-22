function PDE = ePDEsetup1(meshsz,sigma2,l)
if nargin==0
    meshsz=11;
    sigma2=1; l=0.2;
elseif nargin==1
    sigma2=1; l=0.2;
end

% Generate rectangular mesh for [0,1]*[0,1] domain
h = 1/(meshsz-1);
[node,elem] = squaremesh([0,1,0,1],h);

% Compute the coefficients
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
C = sigma2.*exp(-pdist2(center,center).^2/(2*l^2));
[evc,ev] = eigs(C,20);
rtev_KL = sqrt(diag(ev)); evc_KL = evc;

% Boundary conditions
bdFlag = setboundary(node,elem,'Dirichlet','(y==1) | (y==0)');
[bdNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
freeNode = find(~isBdNode);

% output
PDE.node = node; PDE.elem = elem;
PDE.freeNode = freeNode; PDE.bdNode = bdNode;
PDE.rtev_KL = rtev_KL; PDE.evc_KL = evc_KL;
PDE.center = center;
