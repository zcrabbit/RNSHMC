function y = predict(obj,x)
a = exp(-bsxfun(@rdivide,sq_dist(x',obj.center{1}'),obj.width{1}')/2);


y = a*transpose(obj.wts);