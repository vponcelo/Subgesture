function P = createIniPopul(GenomeLength,~,~,X,params)

global BASELINE;
global STATE;

if isempty(STATE)
%     if ~strcmp(params.Baseline,BASELINE{1})
%         seg0ini = seg0(1:2:end);
%         seg0fi = seg0(2:2:end);
%     end
    
    % number of segments for each individual
    if params.probSeg > 0
        nsegs = randi([params.k0 params.N],1,params.population);  % random
    else
        nsegs = params.N*ones(1,params.population);               % N
    end
%     P = zeros(params.population,GenomeLength);
    P = zeros(params.population,params.N*2+1);
    for i = 1:size(P,1)
        %% this run faster for huge N
        idxj = [0 mod(2:size(P,2),2) == 0]; % start seq indices
        if nsegs(i) < params.N
            idxj(1+nsegs(i)*2+1:end) = 0;   % odds indices
            P(i,1+nsegs(i)*2+1:end) = inf;  % odds infinity
        end
        P(i,idxj==1) = randi([1 length(X)-params.nmax+1],1,sum(idxj)); % start seq values
        idxj = P(i,1:sum(idxj)*2+1) == 0;                         % length seq indices
        P(i,idxj) = randi([params.nmin params.nmax],1,sum(idxj)); % length seq
        for j = 2:2:size(P,2)-1
            in = P(i,j);    % choose one non-infinity start
            if in < inf     
                fi = in + P(i,j+1) - 1;
                if ~any(any(X(in:fi,:)))    % check non-zero segment
                    P(i,j:j+1) = P(i,nsegs(i)*2:nsegs(i)*2+1);
                    P(i,nsegs(i)*2:nsegs(i)*2+1) = inf;
                    nsegs(i) = nsegs(i) - 1;
                end
            end
        end
        
        % Sort P
        Peven = P(i,2:2:nsegs(i)*2+1);
        Podd = P(i,3:2:nsegs(i)*2+1);        
        Psort = reshape(sortrows([Peven; Podd]')',1,nsegs(i)*2);
        Psort = [Psort P(i,1+nsegs(i)*2+1:end)];        
        % assign and add k to the sorted population
        P(i,:) = [randi([params.k0 nsegs(i)-1]) Psort]; 
        
        %% Old implementation for different baselines (much slower)
%         in = 1;
%         for j = 1:size(P,2)
%             if mod(j,2) == 0
%                 if nsegs >= round((j-1)/2)
%                     % generate random population
%                     P(i,j) = randi([1 length(X)-params.nmax+1]);
%                     if strcmp(params.Baseline,BASELINE{3})
%                         % generate fixed population 
%                         P(i,j) = in;
%         %             elseif strcmp(params.Baseline,BASELINE{2})
%         %                 % generate derivate population
%         %                 [~,pos] = min(abs(seg0ini-P(i,j)));
%         %                 P(i,j) = seg0ini(pos(1));
%         % %                 P(i,j) = seg0ini(round(j/2));  % comment this line when performing full GA 
%                     end
%                 else
%                     P(i,j) = inf;                
%                 end
%             else
%                 if nsegs >= round((j-1)/2)
%                     % generate random population 
%                     P(i,j) = randi([params.nmin params.nmax]);
%                     if strcmp(params.Baseline,BASELINE{3})
%                         % generate fixed population 
%                         if in+fr_fixed-1 <= size(X,1) && j < size(P,2)
%                             P(i,j) = fr_fixed;
%                         else
%                             P(i,j) = size(X,1)-in+1;
%                         end
%                         in = in+fr_fixed;
%         %             elseif strcmp(params.Baseline,BASELINE{2})
%         %                 % generate derivate population 
%         %                 [~,pos] = min(abs(seg0fi-(P(i,j-1)+P(i,j)-1)));
%         %                 if seg0fi(pos(1)) < P(i,j-1)
%         %                     if pos(1)+1 > length(seg0fi);
%         %                         [~,posIn] = min(abs(seg0ini-P(i,j-1)));
%         %                         fi = seg0ini(posIn);
%         %                     else
%         %                         fi = seg0fi(pos(1)+1);
%         %                     end
%         %                 else
%         %                     fi = seg0fi(pos(1));
%         %                 end
%         %                 P(i,j) = fi-P(i,j-1)+1;
%         % %                 P(i,j) = seg0fi(round(j/2-1))-seg0ini(round(j/2-1))+1;  % comment this line when performing full GA     
%                     end
%                 else
%                     P(i,j) = inf;                
%                 end
%             end
%         end
%         % Sort P
%         [P(i,2:2:nsegs(i)*2+1),is] = sort(P(i,2:2:nsegs(i)*2+1),'ascend');
%         lengths = P(i,3:2:nsegs(i)*2+1);
%         for j = 1:length(is)
%             P(i,2*j+1) = lengths(is(j));
%             if j > 1
%                 if P(i,2*j) == P(i,2*(j-1))
%                     if P(i,2*j+1) < P(i,2*(j-1)+1)
%                         next = P(i,2*(j-1)+1);
%                         P(i,2*(j-1)+1) = P(i,2*j+1);
%                         P(i,2*j+1) = next;
%                     end
%                 end
%             end
%         end
    
    end
    if strcmp(params.mType,'modelSM1') || strcmp(params.mType,'allSM1')
        mnsegs = randi([params.k0 params.N0],1,params.population);
        mk = zeros(1,length(mnsegs));
        for i = 1:length(mk)
            mk(i) = randi([params.k0-1 mnsegs(i)-1]);
        end
        P = [P mnsegs' mk'];
    end
else
    P = STATE.Population;
end
if size(P,2) ~= GenomeLength
    error('createIniPopul:badGenomes','Population genomes are incorrect');
end
