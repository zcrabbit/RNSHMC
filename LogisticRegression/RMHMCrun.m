%% RMHMC 
%% a9a60d
nLeaps = 4;
stepsize = 0.5;
nrep = 5;

for j = 1:nrep
    [betaPosterior,Times,acprat] = BLR_RMHMC(stepsize,nLeaps,'a9a60d');
end
