% Dirichlet Boundary function
function z = g_D(p)
z = p(:,2) + (-1).^p(:,2) .* p(:,1);
%z = cos(2*pi*p(:,1));
end