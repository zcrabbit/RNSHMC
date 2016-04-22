function contourplot(f,axes)
grid = linspace(axes(1),axes(2),100);
[X,Y] = meshgrid(grid,grid);
Z = zeros(size(X));
Z(:) = f([X(:),Y(:)]);
[v, h] = contourf(X,Y,exp(-Z)); %shading flat;
set(h,'Edgecolor','none');
set(gcf, 'colormap', gray);
colormap(flipud(colormap));