function plotPerfMethods(baseline,v1,t1,v2,t2,v3,t3,v4,t4,v5,t5)
    % function that plots the performance for the different methods over 10 runs of the GA

    if ~exist('baseline','var')
        error('plotExecs:errVar','Baseline varaible must exist');
    end
    x = 1:1:10;
    baseline = repmat(baseline,1,10);
    figure,
    subplot(3,2,1);plot(x,baseline,'k-',x,min(v1),'bv',x,max(v1),'b^',x,mean(v1),'b-',x,min(t1),'rv',x,max(t1),'r^',x,mean(t1),'r-')
    hold on
    title('HMM');
    xlabel('Generation');
    ylabel('Score');
    legend('VideoDarwin','maxVal','minVal','meanVal','maxTest','minTest','meanTest','Location','SE');
    subplot(3,2,2);plot(x,baseline,'k-',x,min(v2),'bv',x,max(v2),'b^',x,mean(v2),'b-',x,min(t2),'rv',x,max(t2),'r^',x,mean(t2),'r-')
    title('SRM-HMM');
    xlabel('Generation');
    ylabel('Score');
    legend('VideoDarwin','maxVal','minVal','meanVal','maxTest','minTest','meanTest','Location','SE');
    subplot(3,2,3);plot(x,baseline,'k-',x,min(v3),'bv',x,max(v3),'b^',x,mean(v3),'b-',x,min(t3),'rv',x,max(t3),'r^',x,mean(t3),'r-')
    title('DTW');
    xlabel('Generation');
    ylabel('Score');
    legend('VideoDarwin','maxVal','minVal','meanVal','maxTest','minTest','meanTest','Location','SE');
    subplot(3,2,4);plot(x,baseline,'k-',x,min(v4),'bv',x,max(v4),'b^',x,mean(v4),'b-',x,min(t4),'rv',x,max(t4),'r^',x,mean(t4),'r-')
    title('SRM-DTW');
    xlabel('Generation');
    ylabel('Score');
    legend('VideoDarwin','maxVal','minVal','meanVal','maxTest','minTest','meanTest','Location','SE');
    subplot(3,2,5);plot(x,baseline,'k-',x,min(v5),'bv',x,max(v5),'b^',x,mean(v5),'b-',x,min(t5),'rv',x,max(t5),'r^',x,mean(t5),'r-')
    title('SRM-SVM');
    xlabel('Generation');
    ylabel('Score');
    legend('VideoDarwin','maxVal','minVal','meanVal','maxTest','minTest','meanTest','Location','SE');
    hold off
end