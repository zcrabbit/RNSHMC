function z = rectifier_prime(a) 
%z = a>0;
z = 1./(1+exp(-a));
end