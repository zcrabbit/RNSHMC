function grad = gradient(obj,Xtr) 
n = bsxfun(@plus,Xtr*obj.wts{1}',transpose(obj.biases{1}));
grad = bsxfun(@times,obj.wts{2},rectifier_prime(n))*obj.wts{1};