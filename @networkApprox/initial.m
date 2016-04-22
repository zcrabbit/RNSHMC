function obj = initial(obj,sizes)
obj.numoflayers = length(sizes); obj.sizes = sizes;
obj.biases = cell(obj.numoflayers-1,1);
obj.wts = cell(obj.numoflayers-1,1);
for numlayer = 1:obj.numoflayers-1
    obj.biases{numlayer} = randn(sizes(numlayer+1),1);
    obj.wts{numlayer} = randn(sizes(numlayer+1),sizes(numlayer))/sqrt(sizes(numlayer));
end