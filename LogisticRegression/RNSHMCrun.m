%% RNSHMC 
%% Logistic Regression
nLeaps = [6 10 45];
stepsize = [0.045 0.012 0.012];
s = [2000 2500 1000];
Data = {'simulation' 'a9a60d' 'bankmarket'};
nrep = 5;

for i = 1:length(Data)
    for j = 1:nrep
        [betaPosterior,Times,acprat] = BLR_RNSHMC(s(i),stepsize(i),nLeaps(i),Data{i});
    end
end

    