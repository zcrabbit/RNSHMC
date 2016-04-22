function z = rectifier(a)
%z = a.*(a>0);
z = log(1+exp(a));