function z = f(x)
z = 4*x(:,1).*(pi^2*cos(2*pi*x(:,1)))+2*pi*sin(2*pi*x(:,1));
end
