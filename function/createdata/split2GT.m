function [Xtrain,Xval] = split2GT(X,Y)
% split X into GT segements specified in Y

% Input:
%   X: data sequences
%   Y: sequence labels
% Output:
%   Xtrain: Splitted training sequences 
%   Xval: Splitted validation sequences 

Xaux = cell(1,length(Y));
for j = 1:length(X)
    Xaux{j} = cell(1,length(Y{j}.L));
    fi = 1;
    for i = 1:length(Y{j}.L)
        if i == length(Y{j}.L)
            fi = 0;
        end 
        Xaux{j}{i} = X{j}(Y{j}.seg(i):Y{j}.seg(i+1)-fi,:);    
    end    
end
Xtrain = Xaux{1};
Xval = Xaux{2};
