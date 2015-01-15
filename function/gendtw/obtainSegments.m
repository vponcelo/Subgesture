function [Xtrain,Xval,Y,I] = obtainSegments(X,Y,nSegs,nmin,nmax)
% Obtain initial sequence segments 
%   Output:
%       Y: Learning labels (adding the initial segmentation)
%       Xtrain: Training segmented sequences
%       Xval: Validation segmented sequences
%       I: Individual for the genetic algorithm
%   Input:
%       X: Learning sequences to obtain the segments
%       Y: Learning labels
%       nSegs: number of segments
%       nmin:  minimum subsequence width
%       nmax: maximum sequence width

I = zeros(1,2*nSegs);

Xout = cell(1,length(Y));
for j = 1:length(Y)
    Xout{j} = cell(1,nSegs);
    Y{j}.L0 = 1:nSegs;
    Y{j}.seg0 = []; 
    for i = 1:nSegs
        zeroX = true;
        while zeroX
            in = randi([1 length(X{j})-nmax]);
            fi = in + randi([nmin nmax]) - 1;
            if ~any(any(X{j}(in:fi,:)))
                Xout{j}{i} = X{j}(in:fi,:);
                zeroX = false;
            end                
        end        
        Y{j}.seg0 = [Y{j}.seg0 in fi];
        idx = i*2;
        I(idx-1) = in;
        I(idx) = fi;
    end    
    
    %% segmentación inicial etiquetada
%     i = 1; seg0 = i;
%     while i <= length(X{j}) 
%         if i > length(X{j})-nmax
%             n = length(X{j})-i+1;
%         else
%             n = randi([nmin nmax]);
%         end
%         seg0 = [seg0 i+n-1];
%         i = seg0(end) + 1;    
%     end
%     Y{j}.seg0 = seg0;
%     Y{j}.L0 = crossvalind('Kfold',length(seg0)-1,nSegs)';
end
Xtrain = Xout{1};
Xval = Xout{2};

