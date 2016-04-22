function grad = gradient(obj,x) 
a = exp(-bsxfun(@rdivide,sq_dist(obj.center{1}',x'),obj.width{1})/2);
Dist = bsxfun(@rdivide,bsxfun(@minus,obj.center{1},x),obj.width{1});
grad = obj.wts * bsxfun(@times,Dist,a);