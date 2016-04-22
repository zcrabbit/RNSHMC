function h = forceplot(gf,axes,color)
grid = linspace(axes(1),axes(2),20);
[X,Y] = meshgrid(grid,grid);
grad = -gf([X(:) Y(:)]);
h = quiver(X(:),Y(:),grad(:,1),grad(:,2),color);

