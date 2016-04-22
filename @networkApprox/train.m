function [obj,S,pinS,wts,The] = train(obj,Xtr,Ytr,epsilon)
%if nargin < 4, flag=0;else flag=1; end
if nargin < 4, epsilon = 1e-06; end
%obj.biases{1} = -diag(Xtr*obj.wts{1}');
n = bsxfun(@plus,Xtr*obj.wts{1}',transpose(obj.biases{1}));
S = rectifier(n);
% while det(S) ==0
%     obj = initial(obj,obj.sizes);
%     obj.biases{1} = -diag(Xtr*obj.wts{1}');
%     n = bsxfun(@plus,Xtr*obj.wts{1}',transpose(obj.biases{1}));
%     S = sigmoid(n);
% end
% if flag
%     S = [epsilon*eye(size(S,2));S];
%     Ytr = [zeros(size(S,2),1);Ytr];
% end
%pinS = pinv(S);
% The = pinS*pinS';
%v = pinS*Ytr;
%obj.wts{2} = (pinS*Ytr)';
obj.wts{2} = ((S'*S+epsilon*eye(size(S,2)))\(S'*Ytr))';
%obj.wts{2} = (S\Ytr)';
%obj.biases{2} = v(1);
wts = obj.wts{2}';
%d = graddist(obj,Xtr,Ctr);
% dist = [d];
% 
% while d > tol
%     for i = 1:size(Xtr,1)
%         v = obj.wts{2};v(i) = 0;
%         obj.wts{1}(i,:) = 2/obj.wts{2}(i) *(Ctr(i,:)-(v.*sigmoid_prime(n(i,:)))*obj.wts{1});
%         obj.biases{1}(i) = -Xtr(i,:)*transpose(obj.wts{1}(i,:));
%         n(:,i) = Xtr * transpose(obj.wts{1}(i,:))+obj.biases{1}(i);
%         S = sigmoid(n);
% %         while det(S) == 0
% %             obj = initial(obj,obj.sizes);
% %             obj.biases{1} = -diag(Xtr*obj.wts{1}');
% %             n = bsxfun(@plus,Xtr*obj.wts{1}',transpose(obj.biases{1}));
% %             S = sigmoid(n);
% %         end
%         obj.wts{2} = (pinv(S) * Ytr)';
%         d = graddist(obj,Xtr,Ctr);
%         dist = [dist d];
%         figure(1);plot(dist);
%         if d < tol
%             break;
%         end
%     end
% end

        
    
    
    