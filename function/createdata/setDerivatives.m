function X = setDerivatives(X)
% function that add the derivatives of the input data X, containing the
% x,y,z positions of several skeleton joints

global JOINTS;
global COORDS;

if strcmp(COORDS,'world')
    oc = 3;
elseif strcmp(COORDS,'pixel')
    oc = 2;
end

for i = 1:length(X)
    D = [zeros(1,size(X{i},2)); X{i}(2:end,:)-X{i}(1:end-1,:)];
    D2 = [zeros(1,size(X{i},2)); D(2:end,:)-D(1:end-1,:)];
%     D = [zeros(1,size(X{i}(:,1:length(JOINTS)*oc),2)); X{i}(2:end,1:length(JOINTS)*oc)-X{i}(1:end-1,1:length(JOINTS)*oc)];
%     D2 = [zeros(1,size(X{i}(:,1:length(JOINTS)*oc),2)); D(2:end,1:length(JOINTS)*oc)-D(1:end-1,1:length(JOINTS)*oc)];
    X{i} = [X{i} D D2];
end