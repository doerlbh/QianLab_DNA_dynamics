% by Baihan Lin, August 2016

function unilist = unevenmean(list)
% to output the mean of different depth of mean

x = list(1,:);
y = list(2,:);
uniqx = unique(x);
uniqx = sort(uniqx);
ele = length(uniqx);
uniqy = zeros(1,ele);
for k = 1:ele
	hv = (x==uniqx(k));
    uniqy(k) = dot(hv, y)/sum(hv);
end

unilist = [uniqx; uniqy];

end



