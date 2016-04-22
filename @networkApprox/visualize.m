function visualize(obj,layer)
wts = obj.wts{layer};
biases = obj.biases{layer};
if size(wts,2) == 2
    plot(wts(:,1),wts(:,2),'r.');ax=axis;
    for l = 1:size(wts,1)
        refline([-wts(l,1),-biases(l)]/wts(l,2));
    end
    axis(ax);
else
    plot(wts(1,:),wts(2,:),'r.');
    for l = 1:size(wts,2)
        refline([-wts(1,l),-biases(l)]/wts(2,l));
    end
end
