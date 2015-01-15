function newLabels=displaceSequenceIndexes(labels,percentage),
    newLabels={};
    for i=1:length(labels),
        lab=labels{i};
        newLab=[];
        for j=1:size(lab),
            len=lab(2)-lab(1)+1;
            newLab=[newLab ; round(lab+len*percentage)]
        end
        newLabels=[newLabels ; newLab];
    end
