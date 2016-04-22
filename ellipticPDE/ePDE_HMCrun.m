%% ePDE HMC and RNSHMC run
nrep = 5;
stepsize = 0.16;
nLeap = 10;
s = 1000;
%% HMC
for nr = 1:nrep
    [thetaPosterior,Times,acprat] = ePDE_HMC(20,stepsize,nLeap);
end

%% RNS-HMC
for nr = 1:nrep-1
    [thetaPosterior,Times,acprat] = ePDE_RNSHMC(s,stepsize,nLeap);
end