function twotap_analysis
% purpose: analysing the data from the twotap task just like actual
% behavioral task

% detect pairs of "taps". these are trials
% taps are either defining IPI or ITI

%% find all actions, then sort into ITI or IPI
%close all
%clear all
%% plots of ITI vs trials (intact)
plt = 1;
lrn = 1;
div = 2000;
[results_intact,results_lesioned]=twotap_world(plt,lrn,div);

[tap_idx]=find(results_intact.action == 2);

tap_idx = [tap_idx];
dur_tap = diff(tap_idx); %tells you how much time between each tap, but not whether the tap was IPI or ITI

IPI_start = tap_idx(results_intact.l_state(tap_idx(1:end-1))== 25); %index (timepoint) of when IPI started
IPI_end = tap_idx(find(results_intact.l_state(tap_idx(1:end-1))== 25)+1); %index timepoint of IPI end
IPI=((IPI_end-IPI_start)*100)./1000;
IPI(IPI>1.2) = 1.2;

ITI_start = intersect(IPI_end,tap_idx(results_intact.c_state(tap_idx(1:end-1))== 13));
ITI_end = intersect( tap_idx(find(results_intact.l_state(tap_idx(1:end-2))== 25)+2), tap_idx(find(results_intact.c_state(tap_idx(1:end-1))== 13)+1));

%ITI_start=tap_idx(results_intact.c_state(tap_idx(1:end-1))== 13); %index (timepoint) of when ITI started
%ITI_end = tap_idx(find(results_intact.c_state(tap_idx(1:end-1))== 13)+1); %index timepoint of IPI end
ITI=((ITI_end-ITI_start)*100)./1000;
ITI = ITI;

[a,idx]=sort([ITI_start IPI_start]);
ALL_int=[ITI IPI];

for i = 1:length(ALL_int)
    % basically look at whether the next thing is an IPI or ITI (really
    % dumb way to do things)
    if idx(i)<= length(ITI)
        ITI_final(i) = ALL_int(idx(i));
        IPI_final(i) = NaN;
    else
        ITI_final(i) = NaN;
        IPI_final(i) = ALL_int(idx(i));
    end
end

figure; hold on;
subplot 211
plot(IPI_final,'o')
line([0 length(ALL_int)],[.7 .7],'Color','r')
title('IPI')

subplot 212
plot(ITI_final,'o')
line([0 length(ALL_int)],[1.2 1.2],'Color','r')
title('ITI')
xlabel('trials')
ylabel('duration (s)')

figure
[paramDensity, paramMean, paramCV] = plotDensityLearning(IPI_final*1000');

%% plots of ITI vs trials (lesioned)
clear IPI_start
clear  IPI_end
clear IPI_final
clear ITI_final

clear ITI_start
clear  IPI_end
clear ALL_int

[tap_idx]=find(results_lesioned.action == 2);

tap_idx = [tap_idx];
dur_tap = diff(tap_idx); %tells you how much time between each tap, but not whether the tap was IPI or ITI

IPI_start = tap_idx(results_lesioned.l_state(tap_idx(1:end-1))== 25); %index (timepoint) of when IPI started
IPI_end = tap_idx(find(results_lesioned.l_state(tap_idx(1:end-1))== 25)+1); %index timepoint of IPI end
IPI=((IPI_end-IPI_start)*100)./1000;
IPI(IPI>1.2) = 1.2;

ITI_start = intersect(IPI_end,tap_idx(results_lesioned.c_state(tap_idx(1:end-1))== 13));
ITI_end = intersect( tap_idx(find(results_lesioned.l_state(tap_idx(1:end-2))== 25)+2), tap_idx(find(results_lesioned.c_state(tap_idx(1:end-1))== 13)+1));

%ITI_start=tap_idx(results_lesioned.c_state(tap_idx(1:end-1))== 13); %index (timepoint) of when ITI started
%ITI_end = tap_idx(find(results_lesioned.c_state(tap_idx(1:end-1))== 13)+1); %index timepoint of IPI end
ITI=((ITI_end-ITI_start)*100)./1000;ITI=((ITI_end-ITI_start)*100)./1000;

[a,idx]=sort([ITI_start IPI_start]);
ALL_int=[ITI IPI];

for i = 1:length(ALL_int)
    % basically look at whether the next thing is an IPI or ITI (really
    % dumb way to do things)
    if idx(i)<= length(ITI)
        ITI_final(i) = ALL_int(idx(i));
        IPI_final(i) = NaN;
    else
        ITI_final(i) = NaN;
        IPI_final(i) = ALL_int(idx(i));
    end
end

figure; hold on;
subplot 211
plot(IPI_final,'o')
line([0 length(ALL_int)],[.7 .7],'Color','r')
title('IPI')

subplot 212
plot(ITI_final,'o')
line([0 length(ALL_int)],[1.2 1.2],'Color','r')
title('ITI')
xlabel('trials')
ylabel('duration (s)')

figure

[paramDensity, paramMean, paramCV] = plotDensityLearning(IPI_final*1000');



end