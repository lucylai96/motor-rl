function tapDistances
%close all
% purpose: this function quantifies the Euclidian distance between tap
% trajectories

% a goal is to be able to visualize how taps are explored in a space for% a goal is to be able to visualize how taps are explored in a space for
% lesioned vs nonlesioned animals

%% adding dependencies
addpath('/Volumes/Google Drive/Team Drives/MC Learning Project/Matlab')
addpath('/Volumes/Google Drive/Team Drives/MC Learning Project/Matlab/output/')
addpath('/Volumes/Google Drive/Team Drives/MC Learning Project/Matlab/gerald new code/')

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


%% position AND velocity (taps)

% 16889 trials
%tap_exploration('GP1829_001_061_taps_pv_aah.mat') %naive partial lesion but seemed to learn IPI kind of (after 5000 trials), and ITI kind of!!
%tap_exploration('GP1840_001_059_taps_pv_aah.mat')% learned around 5k trials

% 16869 trials
%tap_exploration('GP1830_001_084_taps_pv_aah.mat')% this guy was switched to auditory feedback but before then did not learn ANY structure

% 20321 trials
%tap_exploration('GP1831_001_045_taps_pv_aah.mat')% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI

% 17084 trials
%tap_exploration('GP1832_001_070_taps_pv_aah.mat') %learned around 13k trials

% 16751 trials
tap_exploration('GP1838_001_099_taps_pv_aah.mat')% learned around 5k trials

% 17800 trials
tap_exploration('GP1840_001_059_taps_pv_aah.mat')% learned around 5k trials

end

function tap_exploration(rat,s)
%% get PC dimensions
analysis = load(rat);
analysis = reconstructTrajectories(analysis);

struct_unzip(analysis);
[dims,plt,idx]=vis(analysis,2);
taptype = taptype(idx);

[dims,plt,idx]=vis(analysis,3); 
tapnum = tapnum(idx);
close

[dims,plt,idx,rewards]=vis(analysis,4);
rewards = rewards(idx);
close

dims = dims(taptype>0,:);% take out 1-tap trials
rewards = rewards(taptype>0);

rewards = rewards(1:end-2);
taptype =taptype(taptype>0);
tapnum =tapnum(taptype>0);
%% what is the distance between tap t and tap t+1
for i = 1:2:length(dims)-2 % for each 1st tap in a 2-tap sequence, how similar/different were they across time
    dists_tap(i) = norm(dims(i+2,:)-dims(i,:));
end

for i = 2:2:length(dims)-1 % for each 2nd tap in a 2-tap sequence, how similar/different were they across time
    dists_tap(i) = norm(dims(i+2,:)-dims(i,:));
end

figure;
subplot 411; hold on;
plot(dists_tap(1:2:1000),'.')
plot(movmean(dists_tap(1:2:1000),10),'LineWidth',1)
prettyplot
subplot 412; hold on;
plot(dists_tap(2:2:1000),'.')
plot(movmean(dists_tap(2:2:1000),10),'LineWidth',1)
prettyplot


subplot 413; hold on;
plot(dists_tap(1:2:end),'.')
plot(movmean(dists_tap(1:2:end),50),'LineWidth',1)
prettyplot
subplot 414; hold on;
plot(dists_tap(2:2:end),'.')
plot(movmean(dists_tap(2:2:end),50),'LineWidth',1)
prettyplot

%% how does reward influence the tap distances at beginning (tr = # trials to look at)
tr = 2000;
% figure;
% plot(rewards(1:tr),dists_tap(1:tr),'ko')
% [c,p]=corrcoef(rewards(1:tr),dists_tap(1:tr))

figure; hold on; subplot 131
boxplot(dists_tap(1:tr),rewards(1:tr))

xlabel('rewards')
ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
title(strcat('First ',num2str(tr),' trials'))
prettyplot
axis square

%% how does reward influence the tap distances at end (tr = # trials to look at)
tr2 = 2000;
% figure;
% scatter(rewards(end-tr2:end),dists_tap(end-tr2:end))
% [c,p]=corrcoef(rewards(end-tr2:end),dists_tap(end-tr2:end))

subplot 132; hold on;
boxplot(dists_tap(end-tr2:end),rewards(end-tr2:end))

%xlabel('rewards')
%ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
title(strcat('Last ',num2str(tr2),' trials'))
prettyplot
axis square

%% boxplots
subplot 133; hold on;
boxplot(dists_tap,rewards)
%text(5, max(cts.Values)+.1*max(cts.Values),strcat('mean: ',num2str(mean(dists_tap(rewards==i)))))
%text(5, max(cts.Values),strcat('std:',num2str(std(dists_tap(rewards==i)))))
%title(strcat('R=',num2str(i), ' | mean: ',num2str(mean(dists_tap(rewards==i))),' | std: ',num2str(std(dists_tap(rewards==i)))))
title('All trials')
%xlabel('rewards')
%ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
prettyplot
axis square
suptitle(rat(1:6))

set(gcf,'Position',[100 450 600 250]);

%% ONLY FIRST TAPS

%% how does reward influence the tap distances at beginning (tr = # trials to look at)
tr = 4000;
% figure;
% plot(rewards(1:tr),dists_tap(1:tr),'ko')
% [c,p]=corrcoef(rewards(1:tr),dists_tap(1:tr))

figure; hold on; subplot 131
boxplot(dists_tap(1:2:tr),rewards(1:2:tr))

xlabel('rewards')
ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
title(strcat('First ',num2str(tr/2),' trials'))
prettyplot
axis square

%% how does reward influence the tap distances at end (tr = # trials to look at)
tr2 = 4000;
% figure;
% scatter(rewards(end-tr2:end),dists_tap(end-tr2:end))
% [c,p]=corrcoef(rewards(end-tr2:end),dists_tap(end-tr2:end))

subplot 132; hold on;
boxplot(dists_tap(end-tr2-1:2:end-1),rewards(end-tr2-1:2:end-1))

%xlabel('rewards')
%ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
title(strcat('Last ',num2str(tr2/2),' trials'))
prettyplot
axis square

%% boxplots
subplot 133; hold on;
boxplot(dists_tap(1:2:end),rewards(1:2:end))
%text(5, max(cts.Values)+.1*max(cts.Values),strcat('mean: ',num2str(mean(dists_tap(rewards==i)))))
%text(5, max(cts.Values),strcat('std:',num2str(std(dists_tap(rewards==i)))))
%title(strcat('R=',num2str(i), ' | mean: ',num2str(mean(dists_tap(rewards==i))),' | std: ',num2str(std(dists_tap(rewards==i)))))
title('All trials')
%xlabel('rewards')
%ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
prettyplot
axis square
suptitle(rat(1:6))

set(gcf,'Position',[100 150 600 250]);


%% ONLY SECOND TAPS

%% how does reward influence the tap distances at beginning (tr = # trials to look at)
tr = 4000;
% figure;
% plot(rewards(1:tr),dists_tap(1:tr),'ko')
% [c,p]=corrcoef(rewards(1:tr),dists_tap(1:tr))

figure; hold on; subplot 131
boxplot(dists_tap(2:2:tr),rewards(2:2:tr))

xlabel('rewards')
ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
title(strcat('First ',num2str(tr/2),' trials'))
prettyplot
axis square

%% how does reward influence the tap distances at end (tr = # trials to look at)
tr2 = 4000;
% figure;
% scatter(rewards(end-tr2:end),dists_tap(end-tr2:end))
% [c,p]=corrcoef(rewards(end-tr2:end),dists_tap(end-tr2:end))

subplot 132; hold on;
boxplot(dists_tap(end-tr2:2:end),rewards(end-tr2:2:end))

%xlabel('rewards')
%ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
title(strcat('Last ',num2str(tr2/2),' trials'))
prettyplot
axis square

%% boxplots
subplot 133; hold on;
boxplot(dists_tap(2:2:end),rewards(2:2:end))
%text(5, max(cts.Values)+.1*max(cts.Values),strcat('mean: ',num2str(mean(dists_tap(rewards==i)))))
%text(5, max(cts.Values),strcat('std:',num2str(std(dists_tap(rewards==i)))))
%title(strcat('R=',num2str(i), ' | mean: ',num2str(mean(dists_tap(rewards==i))),' | std: ',num2str(std(dists_tap(rewards==i)))))
title('All trials')
%xlabel('rewards')
%ylabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
prettyplot
axis square
suptitle(rat(1:6))

set(gcf,'Position',[500 450 600 250]);

%% distribution of unrewarded transitions and rewarded transitions
figure; hold on;
histogram(dists_tap(rewards<1),'FaceColor','r')
histogram(dists_tap(rewards>1),'FaceColor','g')
legend('unrewarded','rewarded')
prettyplot

xlabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
ylabel('# transitions')
suptitle(rat(1:6))

%std(dists_tap(rewards<1)) %unrewarded
%std(dists_tap(rewards>1)) %rewarded

%mean(dists_tap(rewards<1)) %unrewarded
%mean(dists_tap(rewards>1)) %rewarded


% figure; hold on;
% for i = 0:5
%     subplot (2,3,i+1)
%     cts =histogram(dists_tap(rewards==i),50)
%     %text(5, max(cts.Values)+.1*max(cts.Values),strcat('mean: ',num2str(mean(dists_tap(rewards==i)))))
%     %text(5, max(cts.Values),strcat('std:',num2str(std(dists_tap(rewards==i)))))
% title(strcat('R=',num2str(i), ' | mean: ',num2str(mean(dists_tap(rewards==i))),' | std: ',num2str(std(dists_tap(rewards==i)))))
%
%
% end
% equalabscissa(2,3)
% xlabel('|| tap_{1,2}(t) - tap_{1,2}(t+1) || ')
% ylabel('# transitions')
%


%% take any given tap #, how far or close is it from taps that are temporally far or close in distance?
%only first taps
for dis = 2:2:3000 %50 trials on either side
    clear dists_tap
    
     dists_tap =vecnorm(dims(1+dis:2:length(dims),:)-dims(1:2:length(dims)-dis,:),2,2)';
%     for i = 1:2:length(dims)-dis % for each 1st tap in a 2-tap sequence, how similar/different were they across time
%         dists_tap(i) = norm(dims(i+dis,:)-dims(i,:));
%     end
%     
%     for i = 2:2:length(dims)-dis+1 % for each 2nd tap in a 2-tap sequence, how similar/different were they across time
%         dists_tap(i) = norm(dims(i+dis,:)-dims(i,:));
%     end
   dists_tap = dists_tap(1:10000);
    % at the end of this loop, dists_tap is a vector of distances of taps
    % that are (dis) apart
    mega_dists_tap(dis,:) = dists_tap;
end
 mega_dists_tap =  mega_dists_tap(2:2:end,:);

 % each column is a tap
 %each row is the distance it is from the tap(r) away from it
 
 figure;subplot 211;plot(mega_dists_tap(:,1:10)) %any relationship between distance and how far
 subplot 212;plot(mega_dists_tap(:,end-10:end)) %any relationship between distance and how far
xlabel('dis')
ylabel('|| tap(t) - tap(t+dis) || ')
equalabscissa(2,1)

 figure;subplot 211;plot(mega_dists_tap(:,110)) %any relationship between distance and how far
 subplot 212;plot(mega_dists_tap(:,9100)) %any relationship between distance and how far
xlabel('dis')
ylabel('|| tap(t) - tap(t+dis) || ')
equalabscissa(2,1)

%figure;plot(mega_dists_tap(:,100))
%
% tap_rewards = analysis.rewards(analysis.tap_trial(:,1)); %which taps were rewarded by how much
% cmap=colormap(flipud(hot(7)));
% cmap=[1 1 1; 0 0 0];
% tap_rewards(tap_rewards>0) = 1;
% tap_rewards (tap_rewards==0) = 0;
%
% for i = 1:7
%     norms(i,:) =  vecnorm(displacement(:,:,i)');
%
% end
%
%
% figure(1); subplot(2,3,s);hold on;axis([1000 5000 1000 5000])
% xlabel('paw1 ||tap||')
% ylabel('paw2 ||tap||')
% x = 100;
% for tr = 1:x:length(norms)-x
%
%     scatter(norms(4,tr:tr+x-1),norms(5,tr:tr+x-1),[],cmap(tap_rewards(tr:tr+x-1)+1,:),'filled'); %paw 1 and paw 2
%     pause(0.001);
% end
%
%
%


end

function exploration(rat)

analysis = load(rat);
analysis=reconstructTrajectories(analysis);

%% plot distance against time

%% calculate distances between points in PCA spaces
%take x and y coordinates and compute displacement

for tr = 1:x:length(norms)-x
    plot(norms(4,tr:tr+x-1),norms(5,tr:tr+x-1),'ro'); %paw 1 and paw 2
    pause(0.001);
end

%just use x coordinates as displacement
%displacement = squeeze(analysis.tap_filt_rc(:,:,1,:));

%just use y coordinates as displacement
%displacement = squeeze(analysis.tap_filt_rc(:,:,2,:));


holder = NaN*ones(1,61,7);
displacement_shift = cat(1, holder,displacement); %(t-1)
displacement_shift2 = cat(1, holder,displacement_shift); %(t-1)

displacement_shift = displacement_shift(1:end-1,:,:);
displacement_shift2 = displacement_shift2(1:end-2,:,:);


differ = displacement_shift-displacement; %difference in one trial ahead and one trial back [(t-1) - t]
differ2 = displacement_shift2-displacement; %difference in one trial ahead and one trial back [(t-1) - t]

figure; hold on;
for i = 1:size(differ,3)
    euclid_dist(i,:) = vecnorm(differ(:,:,i)');
    euclid_dist2(i,:) = vecnorm(differ2(:,:,i)');
    subplot(4,2,i);plot(euclid_dist(i,:),'.');
    
end
subplot(4,2,1); title('ear')
xlabel('taps')
ylabel('distance between tap(t) and tap(t+1)')
subplot(4,2,2); title('elbow')
subplot(4,2,3); title('nose')
subplot(4,2,4); title('paw')
subplot(4,2,5); title('paw2')
subplot(4,2,6); title('elbow2')
subplot(4,2,7); title('lever')



tap_rewards = analysis.rewards(analysis.tap_trial(:,1)); %which taps were rewarded by how much
figure; hold on;
for i = 1:size(differ,3)
    subplot(4,2,i);scatter(tap_rewards, euclid_dist(i,:));
    
end

subplot(4,2,1); title('ear')
xlabel('reward')
ylabel('distance between tap(t) and tap(t+1)')
subplot(4,2,2); title('elbow')
subplot(4,2,3); title('nose')
subplot(4,2,4); title('paw')
subplot(4,2,5); title('paw2')
subplot(4,2,6); title('elbow2')
subplot(4,2,7); title('lever')


%% are taps following a rewarded tap on average more similar than taps
% following unrewarded taps?
rewarded_taps = euclid_dist(tap_rewards>0);
unrewarded_taps = euclid_dist(tap_rewards==0);
nanmean(euclid_dist(tap_rewards>0))
nanmean(euclid_dist(tap_rewards==0))

figure;hold on
histogram(unrewarded_taps)
histogram(rewarded_taps)
legend('unrewarded','rewarded')

ylabel('distance between tap(t) and tap(t+1)')
xlabel('# taps')

%% same analysis but with magnitude of taps

%color taps by magnitude of reward

figure;hold on
histogram(unrewarded_taps)
histogram(euclid_dist(tap_rewards==1))
histogram(euclid_dist(tap_rewards==2))
histogram(euclid_dist(tap_rewards==3))
histogram(euclid_dist(tap_rewards==4))
histogram(euclid_dist(tap_rewards==5))
histogram(euclid_dist(tap_rewards==6))
legend('unrewarded', 'rew=1','rew=2','rew=3','rew=4','rew=5','rew=6')

ylabel('distance between tap(t) and tap(t+1)')
xlabel('# taps')

nanmean(euclid_dist(tap_rewards==0))
nanmean(euclid_dist(tap_rewards==1))
nanmean(euclid_dist(tap_rewards==2))
nanmean(euclid_dist(tap_rewards==3))
nanmean(euclid_dist(tap_rewards==4))
nanmean(euclid_dist(tap_rewards==5))


%% scatter the distance between t and t+1 and t and t+2 and color by reward
figure; hold on;
for k = 0:6
    subplot(4,2,k+1);scatter(euclid_dist(1,tap_rewards==k), euclid_dist2(1,tap_rewards==k))
    
end

subplot(4,2,1); title('unrewarded')
ylabel('distance between tap(t) and tap(t+2)')
xlabel('distance between tap(t) and tap(t+1)')
subplot(4,2,2); title('rew = 1')
subplot(4,2,3); title('rew = 2')
subplot(4,2,4); title('rew = 3')
subplot(4,2,5); title('rew = 4')
subplot(4,2,6); title('rew = 5')
subplot(4,2,7); title('rew = 6')
equalabscissa(4,2)


idx = find(tap_rewards==5);
figure; hold on;
plot(idx,euclid_dist(1,idx),'.');
xlabel('taps')
ylabel('distance between tap(t) and tap(t+1)')
title('rew =5')

figure; hold on;
for i = 1:size(find(tap_rewards>0))
    plot(idx(i),euclid_dist(1,idx(i)),'.');
    pause(.01)
end
plot(idx,euclid_dist(1,idx),'.');



% plot separate, the taps that were rewarded a lot and the ones that were
% not

figure; hold on
for i = 1:size(differ,3)
    subplot(4,2,i);hold on;plot(find(tap_rewards>0),euclid_dist(i,tap_rewards>0),'b.');
    plot(find(tap_rewards==0),euclid_dist(i,tap_rewards==0),'r.');
end


%scatter(vecnorm(displacement(tap_rewards>0,:,1)'), vecnorm(displacement_shift(tap_rewards>0,:,1)'));
%axis square; axis equal;dline;

%figure
%scatter(vecnorm(displacement(tap_rewards==0,:,1)'), vecnorm(displacement_shift(tap_rewards==0,:,1)'));
%axis square; axis equal;dline;

%take only taps that were rewarded a lot (5) and scatter plot the distance
% between the two



hold on; plot(find(rew1==1),euclid_dist(rew1==1),'rx')


end



function [dims, plt,idx,rewards]=vis(analysis,g)

struct_unzip(analysis);

%home = 'G:\OneDrive - Harvard University\RatStructs\OpCon2\tsne';
%cd(home)
%% parameters

dimname = 'pca'; % 'tsne' or 'pca
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


end