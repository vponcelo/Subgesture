function drawAllSkels(X,segs,mov,lhTrayec,rhTrayec)
% draw skeletons drawing the low hand movements (cyan), trajectory
% changes (yellow), or both (black)
% Input:
%   X: input data sequence to be drawn
%   segs: frame segments
%   mov: logical vector where 1 values are the low hand movements
%   lhTrayec: logical vector where 1 values are the left hand changes of
%       the trajectory
%   rhTrayec: logical vector where 1 values are the right hand changes of
%       the trajectory0

global COORDS;

% X = X.*65535;
bgImage = ones(480,640);
Xdraw = reshape([X(:,1:end/2); X(:,end/2+1:end)]',20,size(X,1),2);
dirname = strcat('results/chalearn2013/alignment/',COORDS,'/');
if ~exist(dirname,'dir')
    mkdir(dirname);
end
k = 1;
for i = 1:size(Xdraw,2)
    lHandColor = 'r'; rHandColor = 'r'; 
    if mov(i)
        lHandColor = 'b';
        rHandColor = 'b';
        if lhTrayec(i) 
            lHandColor = 'g';
        end
        if rhTrayec(i);
            rHandColor = 'g';
        end
    else
        if lhTrayec(i) 
            lHandColor = 'y';
        end
        if rhTrayec(i);
            rHandColor = 'y';
        end
    end
    filename = sprintf('%s%d.png',dirname,i);    
    h = showSkeleton(reshape(Xdraw(:,i,:),20,2),bgImage,1,lHandColor,rHandColor);
    if i >= segs(k) && i <= segs(k+1)
        dirname2 = strcat(dirname,sprintf('/seg%d/',k));
        if ~exist(dirname2,'file')
             mkdir(dirname2);
        end
        filename = sprintf('%s%d',dirname2,i);
        if i == segs(k+1)
            k = k + 1;
        end
    end
    if ~exist(filename,'file')
        saveas(h,filename,'png');
    end
    close(h);
end
