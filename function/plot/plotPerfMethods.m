function plotPerfMethods(baseline,v1,t1,v2,t2,v3,t3,v4,t4,v5,t5)

    % function that plots the performance for the different methods over 10 runs of the GA

    if ~exist('baseline','var')
        error('plotExecs:errVar','Baseline variable must exist');
    end
    x = 1:1:10;
    baseline = repmat(baseline,1,10);
    figure,
    subplot(5,2,1); errorbar(x,mean(v1),std(v1),'b'); %boxplot(v1);
    hold on;
    title('Error bar HMM validation');
    hold off;
    subplot(5,2,2); errorbar(x,mean(t1),std(t1),'r'); %boxplot(t1);
    hold on;
    title('Error bar HMM test');
    hold off;
    subplot(5,2,3); errorbar(x,mean(v2),std(v2),'b'); %boxplot(v2);
    hold on;
    title('Error bar SRM-HMM validation');
    hold off;
    subplot(5,2,4); errorbar(x,mean(t2),std(t2),'r'); %boxplot(t2);
    hold on;
    title('Error bar SRM-HMM test');
    hold off;
    subplot(5,2,5); errorbar(x,mean(v3),std(v3),'b'); %boxplot(v3);
    hold on;
    title('Error bar DTW validation');
    ylabel('Mean score and standard deviations');
    hold off;
    subplot(5,2,6); errorbar(x,mean(t3),std(t3),'r'); %boxplot(t3);
    hold on;
    title('Error bar DTW test');
    ylabel('Mean score and standard deviations');
    hold off;
    subplot(5,2,7); errorbar(x,mean(v4),std(v4),'b'); %boxplot(v4);
    hold on;
    title('Error bar SRM-DTW validation');
    hold off;
    subplot(5,2,8); errorbar(x,mean(t4),std(t4),'r'); %boxplot(t4);
    hold on;
    title('Error bar SRM-DTW test');
    hold off;
    subplot(5,2,9); errorbar(x,mean(v5),std(v5),'b'); %boxplot(v5);
    hold on;
    title('Error bar SRM-SVM validation');
    xlabel('Generation');
    hold off;
    subplot(5,2,10); errorbar(x,mean(t5),std(t5),'r'); %boxplot(t5);
    hold on;
    title('Error bar SRM-SVM test');
    xlabel('Generation');
    hold off;
end