function z = rectifier_prime3(a)
f = 1./(1+exp(a));
z = (1-2*f).*f.*(f-1);