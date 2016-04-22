function [u,pu,p2u] = ePDEsolver1(theta,PDE,opt)
if nargin<3
    opt = 0; 
end

% Compute the coefficients
node = PDE.node; elem = PDE.elem;
rtev_KL = PDE.rtev_KL; evc_KL = PDE.evc_KL;
c=-exp(sum(repmat((theta.*rtev_KL)',[size(evc_KL,1),1]).*evc_KL,2));

N = size(node,1); D = length(theta);

% Compute the stiffness matrix and the right hand side
A = assembling(node,elem,c);
%b = rhandside(node,elem,area,@f);

% Boundary conditions
freeNode = PDE.freeNode; bdNode = PDE.bdNode;

u = zeros(size(node,1),1);
u(bdNode) = g_D(node(bdNode,:));
b_D = - A*u;

% Solve the interior node
u(freeNode) = A(freeNode,freeNode)\b_D(freeNode);
%showsolution(node,elem,u);

if opt>=1
    % compute the derivative of coefficients
    pc = (c*rtev_KL') .* evc_KL;
    pu = zeros(N,D);
    pA = cell(1,D);
    for i = 1:D
        % assembling the stiffness matrix for pu
        pA{i} = assembling(node,elem, pc(:,i));
        pb_D = - pA{i}*u;
        
        % Solve the interior node for pu
        pu(freeNode,i) = A(freeNode,freeNode)\pb_D(freeNode);
    end
    if opt==2
        % compute the second derivative of coefficients
        p2c = bsxfun(@times, pc,reshape(bsxfun(@times,evc_KL,rtev_KL'),[],1,D));
        p2u = zeros(N,D,D);
        for i = 1:D
            for j = 1:i
                p2A = assembling(node,elem,p2c(:,i,j));
                pb_D = -pA{i}*pu(:,j)-pA{j}*pu(:,i)-p2A*u;
                p2u(freeNode,i,j) = A(freeNode,freeNode)\pb_D(freeNode);
                p2u(:,j,i) = p2u(:,i,j);
            end
        end
    end
end

% % verifying
% theta_new = theta + 10^(-3)*[0;0;0;1;0;0];
% c_new = -exp(sum(repmat((theta_new.*rtev_KL)',[size(evc_KL,1),1]).*evc_KL,2));
% 
% A_new = assembling(node,elem,c_new);
% u_new = zeros(size(node,1),1);
% u_new(bdNode) = g_D(node(bdNode,:));
% b_D_new = - A_new*u_new;
% 
% % compare pu
% u_new(freeNode) = A_new(freeNode,freeNode)\b_D_new(freeNode);
% figure;showsolution(node,elem,(u_new-u)/10^(-3));





    
    
