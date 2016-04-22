function H = Hessian(obj,x,der)
if nargin < 3
    der = 0;
end

if der == 0
    V = obj.wts{2}'.*rectifier_prime2(obj.wts{1}*x+obj.biases{1});
    H = obj.wts{1}'*(bsxfun(@times,V,obj.wts{1}));
elseif der == 1
    wts = obj.wts{1};
    for i = 1:length(x)
        V = obj.wts{2}'.*rectifier_prime3(obj.wts{1}*x+obj.biases{1}).*wts(:,i);
        H{i} =  obj.wts{1}'*(bsxfun(@times,V,obj.wts{1}));
    end
end
        
end
