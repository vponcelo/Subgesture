function X = getRandomSubsets(X,k)
% obtain k random subsets of X

minSubsLength = 10^(numel(num2str(length(X)))-3);
I = cell(1,k);
for i = 1:k
    r = 1; c = length(X)/k;
    while r < round(minSubsLength/2)       
        r = length(X)/k;
        while c-r < 1 || c+r > length(X)/k
            c = round(mod(rand(1)*minSubsLength*length(X)/k,length(X)/k));
            r = round(mod(rand(1)*minSubsLength*length(X)/k,length(X)/k)); 
        end
    end
    I{i} = X(c-r:c+r,:);
end