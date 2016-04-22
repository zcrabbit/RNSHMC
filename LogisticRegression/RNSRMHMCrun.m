%% RNS-RMHMC 
%% a9a60d
nLeaps = 4;
stepsize = 0.5;
nrep = 5;
s = 2500;

for j = 1:nrep
    [betaPosterior,Times,acprat] = BLR_RNSRMHMC(s,stepsize,nLeaps,'a9a60d');
end