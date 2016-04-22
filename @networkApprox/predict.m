function y = predict(obj,x)
a = x;
for numlay = 1:obj.numoflayers-2
    a = rectifier(bsxfun(@plus, a*transpose(obj.wts{numlay}), transpose(obj.biases{numlay})));
end

y = a*transpose(obj.wts{numlay+1}) + obj.biases{numlay+1};