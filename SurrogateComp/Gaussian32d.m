%% Compare full GP, additive NN, FITC and RBF NN on a Challenging 32D Gaussian model.

% number of training data
ndata = [300 600 900 1500 2100 3000];
methods = {'GP','NN','FITC','RBF'};
load Data/Sigma;

nrep = 5;

for nr = 1:nrep
    for i = 1:length(methods)
        for nd = ndata
            [thetasaved,Time,acprat,TrainTime] = NarrowGaussian_HMC(Sigma,nd,methods{i});
            CurTime = fix(clock);
            save(['Results/' methods{i} '_' num2str(nd) 'n_' num2str(CurTime) '.mat'], 'thetasaved','Time','acprat','TrainTime');
        end
    end
end
