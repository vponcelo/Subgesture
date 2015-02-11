function P = transform_predictions(predictions)
%   Will transform predictions into csv files to evaluate the performance
%   under the standard protocol of the challence
    P=(reshape(predictions,3,length(predictions)./3))';
end