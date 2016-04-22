%% ePDE HMC run
nrep = 5;
stepsize = 0.16;
nLeap = 10;
s = 1000;
%% HMC
for nr = 1:nrep
    [thetaPosterior,Times,acprat] = ePDE_HMC(stepsize,nLeap);
end
