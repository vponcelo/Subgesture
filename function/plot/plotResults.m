function plotResults(pTrain,pTest,hmmE,s,k)

figure,
if pTrain ~= -1,
    subplot(2,2,1);hist(pTrain);title(sprintf('Histogram of Probabilities for Training Data (%d states)',s));
    xlabel('Probability values');
    ylabel('Number of sequences');
end
if hmmE ~= -1,
    subplot(2,2,2);imagesc(hmmE);title('Labels distribution for each State');
end
initPos = 1;
endPos = initPos+length(pTest)-1;
subplot(2,2,3);hist(pTest(initPos:endPos));title(sprintf('Histogram of Probabilities for Test Data (%d states)',s));
xlabel('Probability values');
ylabel('Number of sequences');

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

mainTitle = sprintf('Results using fold %d',k);
text(0.5, 1, strcat('\bf ',mainTitle),'HorizontalAlignment','center','VerticalAlignment', 'top')

% plotIdx=plotIdx+1;
% initPos=endPos+1;
% endPos=initPos+testLabels{i}(2,2)-1;
% subplot(4,5,plotIdx);hist(pTest{i}(initPos:endPos));title(sprintf('Test Data 10%% displaced (%d states)',s));
% plotIdx=plotIdx+1;
% endPos=initPos+testLabels{i}(3,2)-1;
% subplot(4,5,plotIdx);hist(pTest{i}(initPos:endPos));title(sprintf('Test Data 20%% displaced (%d states)',s));
% plotIdx=plotIdx+1;
% endPos=initPos+testLabels{i}(4,2)-1;
% subplot(4,5,plotIdx);hist(pTest{i}(initPos:endPos));title(sprintf('Test Data 30%% displaced (%d states)',s));
% plotIdx=plotIdx+1;