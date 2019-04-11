function tap_regulation
%close all
% purpose: this function examines how taps are regulated wrt reward

%% adding dependencies
addpath('/Volumes/GoogleDrive/Team Drives/MC Learning Project/Matlab')
addpath('/Volumes/GoogleDrive/Team Drives/MC Learning Project/Matlab/output/')
addpath('/Volumes/GoogleDrive/Team Drives/MC Learning Project/Matlab/gerald new code/')

%% a list of all files:
% these are reconstructed from the dimreduced files
% lesioned
%GP1829_001_061_taps_pv_aah
%GP1830_001_084_taps_pv_aah
%GP1831_001_045_taps_pv_aah

% intact
%GP1832_001_070_taps_pv_aah
%GP1838_001_099_taps_pv_aah
%GP1840_001_059_taps_pv_aah

%% all rats
%tap_exploration('cmulti_taps_p_aah.mat')% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI

%%

tap_exploration('GP1830_001_084_taps_pv_aah.mat')% this guy was switched to auditory feedback but before then did not learn ANY structure
tap_exploration('GP1840_001_059_taps_pv_aah.mat')% learned around 5k trials

%% position AND velocity (taps)
%tap_exploration('GP1840_001_059_taps_pv_aah.mat')% learned around 5k trials
%tap_exploration('GP1830_001_084_taps_pv_aah.mat')% this guy was switched to auditory feedback but before then did not learn ANY structure

% 16889 trials
%tap_exploration('GP1829_001_061_taps_pv_aah.mat') %naive partial lesion but seemed to learn IPI kind of (after 5000 trials), and ITI kind of!!

% 16869 trials
%tap_exploration('GP1830_001_084_taps_pv_aah.mat')% this guy was switched to auditory feedback but before then did not learn ANY structure

% 20321 trials
%tap_exploration('GP1831_001_045_taps_pv_aah.mat')% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI

% 17084 trials
tap_exploration('GP1832_001_070_taps_pv_aah.mat') %learned around 13k trials

% 16751 trials
tap_exploration('GP1838_001_099_taps_pv_aah.mat')% learned around 5k trials

% 17800 trials
%tap_exploration('GP1840_001_059_taps_pv_aah.mat')% learned around 5k trials



%% position (taps)

% 16889 trials
tap_exploration('GP1829_001_061_taps_p_aah.mat') %naive partial lesion but seemed to learn IPI kind of (after 5000 trials), and ITI kind of!!

% 16869 trials
tap_exploration('GP1830_001_084_taps_p_aah.mat')% this guy was switched to auditory feedback but before then did not learn ANY structure

% 20321 trials
tap_exploration('GP1831_001_045_taps_p_aah.mat')% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI

% 17084 trials
tap_exploration('GP1832_001_070_taps_p_aah.mat') %learned around 13k trials

% 16751 trials
tap_exploration('GP1838_001_099_taps_p_aah.mat')% learned around 5k trials

% 17800 trials
tap_exploration('GP1840_001_059_taps_p_aah.mat')% learned around 5k trials



%% velocity (taps)
% 16889 trials
tap_exploration('GP1829_001_061_taps_v_aah.mat') %naive partial lesion but seemed to learn IPI kind of (after 5000 trials), and ITI kind of!!

% 16869 trials
tap_exploration('GP1830_001_084_taps_v_aah.mat')% this guy was switched to auditory feedback but before then did not learn ANY structure

% 20321 trials
tap_exploration('GP1831_001_045_taps_v_aah.mat')% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI

% 17084 trials
tap_exploration('GP1832_001_070_taps_v_aah.mat') %learned around 13k trials

% 16751 trials
tap_exploration('GP1838_001_099_taps_v_aah.mat')% learned around 5k trials

% 17800 trials
tap_exploration('GP1840_001_059_taps_v_aah.mat')% learned around 5k trials


%% multi rat
% multiple rats
tap_exploration('multi_trials_v_aah.mat',3)% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI


end

function tap_exploration(rat)
%close all
mode = 'tsne';
binary = 1; %do on binary rewards or not?
skp = 100;
first =0;
second = 1;
a =2;
%% here add different resolutions in ways of parsing the data
%% plot distance against time
analysis = load(rat);
analysis=reconstructTrajectories(analysis);

%trajPC: ntrials/ntaps x nPCs
% 28928 x 69

%analysis = embedTrajectories(behStruct);



%% first plot the first two PCAs

figure;
subplot 331;
struct_unzip(analysis);
% group numss = {'ratnum (for multi)','1 ipi','2 taptype','3 tapnum','4 rewards','5 rewarded','6 trialnum','7 sessnum','8 daynum'};
% {'tap1','holdtime'} 'rewards50' 'rewarded50'

[dims,plt,idx]=vis(analysis,1,mode); axis square; % group 1 is ipi

[dims,plt,idx]=vis(analysis,2,mode); axis square; % group 2
taptype = taptype(idx);
close

[dims,plt,idx]=vis(analysis,3,mode); axis square; % group 3
tapnum = tapnum(idx);
close

[dims,plt,idx,rewards]=vis(analysis,4,mode); axis square; % group 4
rewards = rewards(idx);
close

%% bin the PCs into 10 equal parts
%taptype = taptype(idx);
%tapnum = tapnum(idx);
%rewards = rewards(idx);
dims = dims(taptype>0,:);% take out 1-tap trials
rewards = rewards(taptype>0);
tapnum =tapnum(taptype>0);
taptype =taptype(taptype>0);

if binary
    rewards(rewards>0)=1;
end
figure;
tapNum = watershed_sandbox(dims,13)';
prettyplot
N = max(tapNum);


%% do taps converge over time?

for i = 1:max(tapNum)
    tap_reward(i) = sum(rewards(tapNum==i),1)./size(rewards(tapNum==i),1);
    tap_freq(i,:) = length(tapNum(tapNum==i))./length(tapNum);
    tap_cumsum(i,:) = cumsum(tapNum==i);
end
% figure;
% subplot 221
% bar(tap_reward)
% xlabel('tap kind')
% ylabel('avg reward rate')
% 
% prettyplot
% 
% subplot 222
% bar(tap_freq);
% xlabel('tap kind')
% ylabel('frequency')
% 
% prettyplot
% 
% subplot (2,2,3:4)
% plot(tap_cumsum','LineWidth',2)
% xlabel('taps (time)')
% ylabel('cumulative sum')
% legend(strsplit(num2str([1:max(tapNum)])))
% prettyplot


%% choose 1st or 2nd taps or all
if first == 1
    tapNum = tapNum(1:2:end);
    rewards = rewards(1:2:end);
elseif second == 1
      tapNum = tapNum(2:2:end);
       rewards = rewards(2:2:end);
end

%% plot the variability of tapNums are being sampled through time
%Y = movvar(tapNum,50);
figure;
subplot 211; hold on;
[cts]=histcounts(tapNum(1:1000));
bar(cts./1000)
%plot(movvar(tapNum(1:1000),10),'LineWidth',1)
%plot(movmean(tapNum(1:1000),10),'LineWidth',1)
title('early (first 1000)')
prettyplot
subplot 212; hold on;
[cts]=histcounts(tapNum(end-1000:end));
bar(cts./1000)
%plot(movvar(tapNum(end-1000:end),10),'LineWidth',1)
%plot(movmean(tapNum(end-1000:end),10),'LineWidth',1)
title('late (last 1000)')
prettyplot
xlabel('tap kind')
ylabel('frequency')
% axxx=subplot(3,1,3); hold on;
% plot(tapNum,'.')
% %plot(movvar(tapNum,100),'LineWidth',1)
% %plot(movmean(tapNum,100),'LineWidth',1)
% title('total')
% linkaxes([ax axx axxx],'y')
% legend('tap kind #','mov variance','mov mean')
% prettyplot

suptitle(rat(1:6))
set(gcf,'Position',[100 100 700 600]);

%% calculate tap frequency
figure; hold on;
norm_tap_freq = zeros(floor(length(tapNum)/skp),max(tapNum));
subplot 121
inc = 1;
for i = 1:length(tapNum)-skp
    [cts]=histcounts(tapNum(i:i+skp), [0.5:1:max(tapNum)+0.5]);
    norm_tap_freq(inc,:) = cts./sum(cts);
    
    inc = inc+1;
end
imagesc(norm_tap_freq)
title('tap frequency')
xlabel('tap modes')
ylabel(strcat('x',num2str(skp),' trials'))
c = colorbar
prettyplot

%% calculate normalized reward rate

subplot 122
inc = 1;

norm_rew_rate = zeros(floor(length(tapNum)/skp),max(tapNum));
for i = 1:length(tapNum)-skp
    for tn = 1:max(tapNum)
        temp_rew = rewards(i:i+skp);
        norm_rew_rate(inc,tn)=sum(temp_rew(tapNum(i:i+skp)==tn))/length(temp_rew(tapNum(i:i+skp)==tn));
    end
    inc = inc+1;
end


imagesc(norm_rew_rate)
title('avg reward values')
colormap(gca,flipud(hot))
c = colorbar
prettyplot


% take out really infrequently occuring taps (only take taps that occured
% >10% of the time
norm_tap_freq(norm_tap_freq<0.05) = NaN;
norm_rew_rate(isnan(norm_tap_freq)) = NaN; %same with associated reward values

%% plot small versions of the above

% figure; hold on;
% subplot 121
% imagesc(norm_tap_freq(end-1000:end,:))
% title('tap frequency')
% xlabel('tap modes')
% ylabel(strcat('x',num2str(skp),' trials'))
% c = colorbar
% prettyplot
% subplot 122
% imagesc(norm_rew_rate(1000:end,:))
% title('avg reward values')
% colormap(gca,flipud(hot))
% c = colorbar
% prettyplot
%% correlate the reward values with tap frequencies in time
% figure
% cor = NaN(floor(size(norm_rew_rate,1)/skp),(1000/skp)+1);
% p_cor = NaN(floor(size(norm_rew_rate,1)/skp),(1000/skp)+1);
% id = 1;
% for r = 1:skp:size(norm_rew_rate,1)-(1000/skp) % for each reward value column, correlate with all of tap frequency rows
%     [cor(id,:) p_cor(id,:)] = corr(norm_rew_rate(r,:)',norm_tap_freq(r:r+(1000/skp),:)','rows','complete'); %only look at correlations out to 1000 trials after
%     id = id+1;
% end
% % each row is a new reward value vector
% %imagesc(cor')
% plot(cor(10:15,:)'); hold on; plot(cor(end-10:end-5,:)');
% xlabel(strcat('lag (x',num2str(skp),') trials'))
% ylabel('correlation')
%% are taps that are rewarded more than avg being upregulated? (early)

total_avg_rew=nanmean(norm_rew_rate,2);
norm_avg_rew=(norm_rew_rate)./total_avg_rew;
figure; hold on
for i = 1:max(tapNum)
subplot (2,ceil(max(tapNum)/2),i)
plot(norm_avg_rew(1:2000,i),'LineWidth',2);hold on;plot(norm_tap_freq(1:2000,i)./norm_tap_freq(1,i),'LineWidth',2);
title(strcat('tap mode ',num2str(i)))
prettyplot
if i ==1
legend('avg reward rate','tap frequency')
xlabel('taps (moving window)')
ylabel('frequency/rate')
end
end
suptitle('early')
%% are taps that are rewarded more than avg being upregulated? (late)

figure; hold on
for i = 1:max(tapNum)
subplot (2,ceil(max(tapNum)/2),i)
plot(norm_avg_rew(end-2000:end,i),'LineWidth',2);hold on;plot(norm_tap_freq(end-2000:end,i)./norm_tap_freq(end-2000,i),'LineWidth',2);
title(strcat('tap mode ',num2str(i)))
prettyplot
if i ==1
legend('avg reward rate','tap frequency')
xlabel('taps (moving window)')
ylabel('frequency/rate')
end
end
suptitle('late')
% [acor,lag] = xcorr(norm_avg_rew(end-2000:end,7)',norm_tap_freq(end-2000:end,7)');
%     [~,I] = max(abs(acor));
%     timeDiff = lag(I);       % sensor 2 leads sensor 1 by 350 samples
%     subplot(311); plot(norm_tap_freq(end-2000:end,7)./norm_tap_freq(end-2000,7)); title('s1');
%     subplot(312); plot(norm_avg_rew(end-2000:end,7)); title('s2');
%     subplot(313); plot(lag,acor);
%     title('Cross-correlation between s1 and s2')

%% correlate the reward values with RATIO of tap frequencies in time (how does tap ratio change as a function of reward)

cor = NaN(floor(size(norm_avg_rew,1)/skp),(1000/skp));
p_cor = NaN(floor(size(norm_avg_rew,1)/skp),(1000/skp));
id = 1;
%cor = NaN(floor(length(tapNum)/skp));
figure; hold on;
for r = 1:skp:size(norm_avg_rew,1)-(1000/skp)
    %  for k = 0:length(norm_tap_freq)-r
    norm_tap_freq_ratio = norm_tap_freq./norm_tap_freq(r,:);
    norm_tap_freq_ratio  = norm_tap_freq_ratio +(rand(size(norm_tap_freq_ratio))*0.1);
    %end
    %% early
    if r == skp+1
        subplot 211
        plot(norm_avg_rew(r,:),norm_tap_freq_ratio(r:r+5,:),'o-')
    end
    %% late
    if r == size(norm_avg_rew,1)-5*skp+1
        subplot 212; hold on;
        plot(norm_avg_rew(r,:),norm_tap_freq_ratio(r:r+5,:),'o-')
    end
    
    if r<length(norm_avg_rew)-(1000/skp)
        
        % for i = 0:(1000/skp) %only calculate out to 1000 trials after that reward %length(norm_tap_freq_ratio)-r
        [cor(id,:) p_cor(id,:)] = corr(norm_avg_rew(r,:)', norm_tap_freq_ratio(r+1:r+(1000/skp),:)','rows','complete');
        id = id+1;
        %plot(norm_avg_rew(1,:),norm_tap_freq_ratio(r:end,:),'o-')
    end
    %  end
    
    %size(norm_tap_freq_ratio)
    %clear norm_tap_freq_ratio
end

figure;hold on
axis([ 1 11 -1 1]); 

plot(cor(1:5,:)','b','LineWidth',1); hold on; plot(cor(end/2:end/2 + 5,:)','g','LineWidth',1); plot(cor(end-5:end,:)','r','LineWidth',1);
xlabel(strcat('lag (x',num2str(skp),') trials'))
ylabel('\rho (reward, \Delta tap freq)')

% figure
% for r = 1:size(norm_avg_rew,1)
%     for i = 1:size(norm_tap_freq,1)
%         cor(r,i) = corr(norm_avg_rew(r,:)', norm_tap_freq(i,:)');
%     end
% end
% % each row is a new reward value vector
% plot(cor(1,:)')
% xlabel('lag')
% ylabel('correlation')


%% make the bins into a matrix (NxN)
% these are the N^2 diff kinds of taps

%tapNum = transformMatrix(bins); %the tap "number"

% have to assign each tap to a bin
% T is N^2

%subplot 339; hold on;
% a = 1;
% 
% T = zeros(N^a);
% R = zeros(N^a);
% %T_plt=zeros(N^2,N^2, mod(size(tapNum),100));
% for i = 1:2:size(tapNum,2)-2 %only look at the tap 1--> tap 2 transitions
%     
%     T(tapNum(i), tapNum(i+1)) =  T(tapNum(i), tapNum(i+1))+1; %fill in transition matrix
%     R(tapNum(i), tapNum(i+1)) = R(tapNum(i), tapNum(i+1))+rewards(i);%rewards for that transition
% end
% T(T<0.02*length(tapNum)/N)=0;
% R(T==0)=0;
% T = T./sum(T(:)); %proportion transition matrix
% R = R./sum(R(:)); %proportion transition matrix
% 
% 
% 
% %% shuffling all taps
% T2 = zeros(N^a);
% R2 = zeros(N^a);
% for p = 1:500
%     scr = randperm(length(tapNum));
%     tapNum2= tapNum(scr);
%     rewards_scr = rewards(scr);
%     %disp(tapNum2(1:5))
%     for i = 1:2:size(tapNum2,2)-1
%         T2(tapNum2(i), tapNum2(i+1),p) =  T2(tapNum2(i), tapNum2(i+1))+1; %fill in transition matrix
%         R2(tapNum2(i), tapNum2(i+1),p) = R2(tapNum2(i), tapNum2(i+1))+rewards_scr(i);%rewards for that transition
%         
%     end
%     % T2(:,:,p) = T2(:,:,p)./sumCts;
% end
% 
% TTT=T./prctile(T2,97.5,3); %most conservative estimate
% TTT(TTT<1) =0; % only the significant ratios (>1)
% 
% %log(TTT)
% %TTT=TT./max(T2,[],3); %most conservative estimate
% figure;
% subplot 131
% imagesc(log(TTT))
% imagesc(T)
% 
% xlabel('tap(t+1)')
% ylabel('tap(t)')
% axis square
% c = colorbar;
% c.Label.String =  'Expected:Shuffled Transitions';
% title('Transition Matrix')
% 
% set(gca,'ydir','normal')
% set(gcf,'color','w');
% 
% prettyplot(12)
% % rewards
% RRR=R./prctile(R2,97.5,3); %most conservative estimate
% RRR(RRR<1) =0;
% 
% subplot 132
% imagesc(log(RRR));
% imagesc(R);
% set(gca,'YDir','normal');
% xlabel('tap(t+1)')
% ylabel('tap(t)')
% axis square
% colormap(gca,flipud(hot));
% c = colorbar;
% c.Label.String = 'Expected:Shuffled Transitions';
% title('Rewards Matrix')
% prettyplot(12)
% 
% % subplot 133
% % scatter(reshape(TTT,[1 N^(a*2)]),reshape(RRR,[1 N^(a*2)]));axis square;
% % %TTT(TTT<0) = NaN;
% % %RRR(RRR<0) = NaN;
% % 
% % [r,p]=corr(reshape(TTT,[1 N^(a*2)])',reshape(RRR,[1 N^(a*2)])','row','complete')
% % text(3,1,strcat('r^2:',num2str(r^2)))
% % text(3,.5,strcat('p:',num2str(p)))
% % title('Correlation')
% % prettyplot(12)
% % suptitle(rat(1:6))
% % set(gcf,'Position',[100 100 900 270]);
% 
% 
% %% separate by first and second tap in a trial
% 
% %% do taps converge over time?
% figure;
% subplot 121
% plot(tapNum,'.')
% xlabel('time')
% ylabel('tap kind')
% axis square
% prettyplot
% subplot 122
% histogram(tapNum);
% xlabel('tap kind')
% ylabel('number of taps')
% axis square
% prettyplot
% 
% %% looking at how rewards relates to tapNum overall
% 
% for i = 1:max(tapNum)
%     tap_reward(i) = size(rewards(tapNum==i),1);
%     tap_cumsum(i,:) = cumsum(tapNum==i);
% end
% figure;
% subplot 121
% bar(tap_reward)
% xlabel('tap kind')
% ylabel('cumulative rewards')
% prettyplot
% 
% subplot 122
% plot(tap_cumsum','LineWidth',2)
% xlabel('taps (time)')
% ylabel('cumulative sum')
% legend('1','2','3','4','5','6','7','8','9')
% prettyplot


%figure; scatter(tapNum,rewards);

%% separate by first and second tap in a trial


%struct_unzip(analysis)
%plot(rewards,'.')
end

function [dims, plt,idx,rewards]=vis(analysis,g,dimname)

struct_unzip(analysis);

%home = 'G:\OneDrive - Harvard University\RatStructs\OpCon2\tsne';
%cd(home)
%% parameters

%dimname = 'tsne'; % 'tsne' or 'pca
showVideo = 4; %1 for paw arcs, 2 whole body arcs, 3 for limbs, 4 for limbs+arcs

plot_parts = [4 5 7]; % parts to plot (Fig2/4). default dom/opp paw, ear

plotOpp = false; % plot opp paw and camera

%% color scatter plot by groups
groups = {'ratnum','ipi','taptype','tapnum','rewards','rewarded','trialnum','sessnum','daynum'};
% 'ratnum','taptype','ipi','rewarded','daynum','pdiff'
% 'rewards','tapnum','rewarded50','rewards50','sday','rday'

% define custom grouping
% groups = {'custom'};
% custom = bi2de([tapnum==2, sessnum>prctile(sessnum,50)]);
% customstr = 'Tap# and Session#';

%% show subset of trials
plt = true(size(ipi));
% plt = sessnum>90;
% plt = mod(daynum,10)==9;
% plt = daynum==9;
% plt = sessnum<prctile(sessnum,25) | sessnum>prctile(sessnum,75);
% plt = rewards>0;% & (sessnum<188 | sessnum>205);
% plt = sessnum<166;
% plt = sessnum>45 & ipi>600 & ipi<800;
% plt = (sessnum>30 & sessnum<40);% & rewards==0;

%% default parameters
useRecon = 1; % reconstruct

pnames = {'ear','dom elbow','nose','dom paw','opp paw','opp elbow','lever'};

% heatmap (Fig2) and line graphs (Fig3)
select3 = true; % click selects 3 nearby trials
plotTapDur = 0; % show dashed lines for tap hold times
plotVel = true; % plot velocity instead of position
plotWarped = true; % plot warped velocity traces
anRaw = false; % plot raw analog instead of velocity

% video stuff
repeats = 1; % video replay
slow_factor = 5; % video speed

% k nearest neighbors
recomputeNN = false; % if true, recompute
k=20;
neighbors=[1 2]; % which neighbors to select

cmap = distinguishable_colors(njoints,{'w','k'});

% figure positions
% if strcmp(getenv('COMPUTERNAME'),'GERALD_LAB')
%     fpos = [5   547   512   420; ...
%         6    43   793   415; ...
%         520   546   278   420;...
%         803    40   873   418;...
%         801   427   441   539;...
%         801   202   553   764];
% else
%     fpos = [5   547   512   420; ...
%         6    43   793   415; ...
%         520   546   278   420;...
%         803    40   873   418;...
%         801   427   441   539;...
%         801   202   553   764];
% end

multiFlag = contains(datatype,'multi');
n=length(plot_parts);
%% load appropriate data

if ~exist('traj','var')
    useRecon=1;
elseif tapFlag && ~exist('tap_filt_rc','var')
    useRecon=0;
elseif ~tapFlag && ~exist('traj_warp_rc','var')
    useRecon=0;
end
if useRecon>0 || tapFlag
    plotWarped = true;
end

% store parameters
vnames = who;
fnames = [fieldnames(analysis); 'analysis'; 'vnames'; 'ans'];
params = struct_zip(vnames(~ismember(vnames,fnames)));
analysis.params = params;

if useRecon>0
    disp('using reconstructed trajectories')
    if tapFlag
        tap_filt = tap_filt_rc;
        tap_vel = tap_vel_rc;
        tap_speed = sqrt(nansum(tap_vel.^2,3));
    else
        traj_warp = traj_warp_rc;
        traj_warp_vel = traj_warp_vel_rc;
        traj_warp_speed = sqrt(nansum(traj_warp_vel.^2,3));
    end
end

if tapFlag
    % replace traj variables with tap variables
    time_warp = time_tap;
    traj_warp = tap_filt;
    traj_warp_vel = tap_vel;
    traj_warp_speed = tap_speed;
    warpFactor = ones(size(tap_filt,1),1);
    
    % re-index trial-based variables with taps
    fnames = fieldnames(analysis);
    ntrials = length(ipi);
    for i = 1:length(fnames)
        if size(analysis.(fnames{i}),1)==ntrials
            eval(sprintf('%s=%s(tap_trial(:,1),:,:,:);',fnames{i},fnames{i}));
        end
    end
    tap1 = tap_trial(:,3);
    tap2 = nan(size(tap1));
    tapnum = tap_trial(:,2);
    
    plt = plt(tap_trial(:,1));
else
    groups(strcmp(groups,'tapnum'))=[];
    groups(strcmp(groups,'taptype'))=[];
end

if ~multiFlag
    groups(strcmp(groups,'ratnum'))=[];
end

% recompute trial number in session
if ~exist('sesstrialnum','var')
    sesstrialnum = nan(size(sessnum));
    if ~multiFlag; rats = 1; end
    for r = 1:length(rats)
        if multiFlag
            behStruct = output{r}.behStruct;
        end
        for s = 1:max(sessnum)
            ix = sessnum==s;
            if multiFlag
                ix = ix & ratnum==r;
            end
            sesstrialnum(ix) = round(trialnum(ix)*length(behStruct(s).rewards));
        end
    end
end

if anRaw
    traj_speed(:,:,1,an) = traj_filt(:,:,1,an);
    traj_warp_speed(:,:,1,an) = traj_warp(:,:,1,an);
end

plt = plt(kp);

%g = 1; % group index
%close all
value = 0;
plotted1 = false; % has heatmap been plotted?
plotted2 = false; % have individual trials been plotted
%while value(1)~=113
%% figure 1: scatter plot of embeddings
% h1=figure(1);
%figure
%    WinOnTop(h1);
%  set(h1,'pos',fpos(1,:))

if ~strcmpi(dimname,'pca')
    dims = trajTSNE;
else
    dims = trajPC;
end
groupby = groups{g};
idx = kp(plt);

selectMode = 2; % 1: select points; 2: select trials;
ix = find(rewards>3); % 0: starting selection

switch lower(groupby)
    case 'ipi'
        
        scatter(dims(plt,1),dims(plt,2),[],ipi(idx),'.'); str = 'IPI (ms)'; cmap1 = parula(1201); cmap1(1,:) = 0; colormap(gca,cmap1);
        title(str)
    case {'sessnum','session','sess'}
        
        scatter(dims(plt,1),dims(plt,2),[],sessnum(idx),'.'); str = 'Session #'; colormap(gca,jet(max(sessnum)));
        title(str)
    case {'tap1','holdtime'}
        
        scatter(dims(plt,1),dims(plt,2),[],tap1(idx),'.'); caxis([0 500]); str = 'Hold time (ms)';
        title(str)
    case 'daynum'
        
        scatter(dims(plt,1),dims(plt,2),[],daynum(idx),'.'); str = 'Day'; colormap(gca,jet(max(daynum)));
        title(str)
    case 'trialnum'
        scatter(dims(plt,1),dims(plt,2),[],trialnum(idx),'.'); str = 'Trial number';
        title(str)
    case 'rewards'
        scatter(dims(plt,1),dims(plt,2),[],rewards(idx),'.'); str = 'Rewards'; colormap(gca,jet(6));
        title(str)
    case 'rewarded'
        scatter(dims(plt,1),dims(plt,2),[],rewards(idx)>0,'.'); str = 'Rewarded';
        title(str)
    case 'tapnum'
        scatter(dims(plt,1),dims(plt,2),[],tapnum(idx),'.'); str = 'Tap #'; colormap(gca,jet)
        title(str)
    case 'taptype'
        colors=[0 0 0; 0 0 1; 0.5 0.5 1; 1 0.5 0.5; 1 0 0];
        scatter(dims(plt,1),dims(plt,2),[],colors(taptype(idx)+1,:),'.'); str = 'Tap Type';
        colormap(gca,colors)
        title(str)
    case 'rewards50'
        scatter(dims(plt,1),dims(plt,2),[],rewards50(idx),'.'); str = 'Avg Rewards';
        title(str)
    case 'rewarded50'
        scatter(dims(plt,1),dims(plt,2),[],rewarded50(idx),'.'); str = 'Avg Reward Rate';
        title(str)
    case 'ratnum'
        scatter(dims(plt,1),dims(plt,2),[],ratnum(idx),'.'); str = 'Rat'; colormap(gca,distinguishable_colors(max(ratnum),'w'))
        title(str)
end


xlabel([dimname ' 1'])
c = colorbar;
c.Label.String = str;
ylabel([dimname ' 2'])



%end
end