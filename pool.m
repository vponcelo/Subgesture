function pool(nworkers)
% input:
%   nworkers: maximum number of workers you wish to use.

while ~matlabpool('size') && nworkers > 0
    try
        matlabpool('open',nworkers);
    catch e
        display(e.cause{1}.message);
        nworkers = nworkers - 1;
    end
end