%% RNSHMC 
%% ePDE
nLeaps = 10;
stepsize = 0.16;
s = 1000;
nrep = 5;

for j = 1:nrep
    [thetaPosterior,Times,acprat] = ePDE_RNSHMC(s,stepsize,nLeaps);
end
