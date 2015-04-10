function detSw = getActivations(pred, gt, segs, minOverlap)
% get predicted activations from GT segments

detSw = zeros(1,length(segs)-1);      % logical vector of the label length
d = diff(pred); idxL = find(d~=0)+1;    % get segments of transitions (derivatives)

if ~isempty(idxL)
    if all(pred), idxL = [1 length(d)]; end
    if d(idxL(1)-1) < 0, idxL = [1 idxL]; end
    idxL(2:2:end) = idxL(2:2:end)-1;    % fix the ending position of the negative transitions to avoid getting additional GT segments
    if d(idxL(end)-1) > 0 && mod(length(idxL),2), idxL = [idxL length(pred)]; end
    if mod(length(idxL),2) > 0,
        error('g:transErr','Transition indexs must be pairs of activation segments'); 
    end
    for l = 1:2:length(idxL)        % treat activations
        % get corresponding GT labels from predicted activation segments, if any
        if idxL(l) == segs(end), sLabel = sum(segs < idxL(l)); else sLabel = sum(segs <= idxL(l)); end
        if idxL(l+1) == segs(end), eLabel = sum(segs < idxL(l+1)); else eLabel = sum(segs <= idxL(l+1)); end
        
        if sLabel > 0
            if segs(eLabel+1) == segs(end), f = 0; else f = -1; end    % avoid exceed the ending index of the last segment
            GTfr = gt(segs(sLabel):segs(eLabel+1)-f);   % get gt segments from predicted activation segments
            Predfr = zeros(1,length(GTfr));             % predicted activations from GT segments    
            gaps = eLabel-sLabel; iDetL = [];
            if eLabel == sLabel     % predicted activation segment belongs to the same GT segment
                Predfr(idxL(l)-segs(sLabel)+1:idxL(l+1)-segs(sLabel)+1) = 1;
                iDetL = [iDetL sLabel];
            else                    % predicted activation segment belongs to several GT segments
                iDetL = [iDetL sLabel];
                Predfr(idxL(l)-segs(sLabel)+1:segs(sLabel+1)-segs(sLabel)) = 1;
                for m = 2:gaps      % worst case: more than one GT segments with activations
                    eTmp = sLabel+m;
                    if all(gt(segs(eTmp):segs(eTmp+1)-1))   % assign intra-activations
                        iDetL = [iDetL eTmp];
                        if eLabel == eTmp
                            Predfr(segs(eLabel)-segs(sLabel)+1:end) = 1;
                        else
                            Predfr(segs(eTmp)-segs(sLabel)+1:segs(eTmp+1)-segs(sLabel)+1) = 1;
                        end                            
                    end
                end
            end
            
            % assign predicted (activated) labels
            if sum(GTfr & Predfr)/sum(GTfr | Predfr) > minOverlap
                detSw(iDetL) = 1;
            end
        end
    end
end