function out = visualizeEmbeddings(analysis)
if nargin<1 || isempty(analysis)
    analysis = evalin('caller','analysis');
end
if analysis.compressAndSave
    if ~any(ismember(fieldnames(analysis),{'tap_filt_rc','traj_warp_rc'}))
        analysis = reconstructTrajectories(analysis);
    end
end
struct_unzip(analysis);

%home = 'G:\OneDrive - Harvard University\RatStructs\OpCon2\tsne';
%cd(home)
%% parameters

dimname = 'pca'; % 'tsne' or 'pca
showVideo = 4; %1 for paw arcs, 2 whole body arcs, 3 for limbs, 4 for limbs+arcs

plot_parts = [4 5 7]; % parts to plot (Fig2/4). default dom/opp paw, ear

plotOpp = false; % plot opp paw and camera

%% color scatter plot by groups
groups = {'ratnum','taptype','ipi','rewarded','sessnum','daynum','pdiff','sday'};
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
if strcmp(getenv('COMPUTERNAME'),'GERALD_LAB')
    fpos = [5   547   512   420; ...
            6    43   793   415; ...
          520   546   278   420;...
          803    40   873   418;...
          801   427   441   539;...
          801   202   553   764];
else
    fpos = [5   547   512   420; ...
            6    43   793   415; ...
          520   546   278   420;...
          803    40   873   418;...
          801   427   441   539;...
          801   202   553   764];
end

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

selectMode = 1; % 1: select points; 2: select trials; 
ix = find(rewards>3); % 0: starting selection

%% find nearest neighbor day/session
if recomputeNN || ~exist('nnidx','var') || any(nnidx(~isnan(nnidx))>length(ipi))
    disp('recomputing neighbors...')
    
    if tapFlag;  mu = tap_mu_vel;
    else;  mu = traj_mu_vel;     end
    
    temp = traj_warp_vel; temp = temp-mu; temp(:,:,:,7)=[];
    temp = reshape(temp,size(temp,1),[]);
    NS = createns(temp);
    [temp, temp2] = knnsearch(NS,temp,'K',k+1);    
    temp(:,1) = []; temp2(:,1)=[]; n = size(ipi,1); 
    nnidx = temp; nndist= temp2;
        
%     NS = createns(trajPC);
%     [temp, temp2] = knnsearch(NS,trajPC,'K',k+1);
%     temp(:,1) = []; temp2(:,1)=[]; n = size(ipi,1);    
%     nnidx = nan(n,k); nndist = nan(n,k);
%     nnidx(kp,:) = kp(temp); nndist(kp,:) = temp2;
end

mthresh = nanmedian(nndist(:,1)); % median of 1st neighbor
fprintf('Median distance: %1.1f\n',mthresh)

weights = 1./(1+exp(0.025*(nndist-200)));


pdays = nan(size(nnidx));
pdays(kp,:) = daynum(nnidx(kp,:));
pdays(nndist>mthresh*3)=NaN;

pday = nanmedian(pdays,2);
rday = nanmax(pdays,[],2)-nanmin(pdays,[],2); % range
sday = nanstd(pdays,[],2); % std
pdiff = pday-daynum;

pdiff = sum(pdiff.*weights,2)./sum(weights,2);

if recomputeNN && multiFlag
    groups(strcmp(groups,'pday'))=[];
    groups(strcmp(groups,'pdiff'))=[];
    groups(strcmp(groups,'rday'))=[];
    groups(strcmp(groups,'sday'))=[];
end

%% swap video side
if plotOpp
dom = opp;
opp = setxor(dom,{'L','R'}); opp = opp{1};
paw = 5; elbow=6; paw2=4; elbow2=2;
plot_parts = [5 4 7];

params.dom = dom; 
params.opp = opp;
params.paw = paw;
params.paw2 = paw2;
params.elbow = elbow;
params.elbow2 = elbow2;
params.plot_parts = plot_parts;
analysis.params = params;
end

%% plot tSNE
plt = plt(kp);

g = 1; % group index
close all
value = 0; 
plotted1 = false; % has heatmap been plotted?
plotted2 = false; % have individual trials been plotted
while value(1)~=113
    %% figure 1: scatter plot of embeddings
    h1=figure(1);
    figure(1)
    WinOnTop(h1);
    set(h1,'pos',fpos(1,:))
    
    if ~strcmpi(dimname,'pca')
        dims = trajTSNE;
    else
        dims = trajPC;
    end
    groupby = groups{g};
    idx = kp(plt);
    
    switch lower(groupby)
        case 'ipi'                        
            scatter(dims(plt,1),dims(plt,2),[],ipi(idx),'.'); str = 'IPI (ms)'; cmap1 = parula(1201); cmap1(1,:) = 0; colormap(cmap1);
        case {'sessnum','session','sess'}
            scatter(dims(plt,1),dims(plt,2),[],sessnum(idx),'.'); str = 'Session #'; colormap(jet(max(sessnum)));
        case {'tap1','holdtime'}
            scatter(dims(plt,1),dims(plt,2),[],tap1(idx),'.'); caxis([0 500]); str = 'Hold time (ms)';
        case 'daynum'
            scatter(dims(plt,1),dims(plt,2),[],daynum(idx),'.'); str = 'Day'; colormap(jet(max(daynum)));
        case 'trialnum'
            scatter(dims(plt,1),dims(plt,2),[],trialnum(idx),'.'); str = 'Trial number';
        case 'rewards'
            scatter(dims(plt,1),dims(plt,2),[],rewards(idx),'.'); str = 'Rewards'; colormap(jet(6));
        case 'rewarded'
            scatter(dims(plt,1),dims(plt,2),[],rewards(idx)>0,'.'); str = 'Rewarded';
        case 'tapnum'
            scatter(dims(plt,1),dims(plt,2),[],tapnum(idx),'.'); str = 'Tap #'; colormap(jet)
        case 'taptype'
            colors=[0 0 0; 0 0 1; 0.5 0.5 1; 1 0.5 0.5; 1 0 0]; colormap(colors)
            scatter(dims(plt,1),dims(plt,2),[],colors(taptype(idx)+1,:),'.'); str = 'Tap Type'; 
        case 'rewards50'
            scatter(dims(plt,1),dims(plt,2),[],rewards50(idx),'.'); str = 'Avg Rewards';
        case 'rewarded50'
            scatter(dims(plt,1),dims(plt,2),[],rewarded50(idx),'.'); str = 'Avg Reward Rate';
        case 'ratnum'
            scatter(dims(plt,1),dims(plt,2),[],ratnum(idx),'.'); str = 'Rat'; colormap(distinguishable_colors(max(ratnum),'w'))
        
        % nearest neighbor -based groupings
        case 'pday'
            scatter(dims(plt,1),dims(plt,2),[],pday(idx),'.'); str = 'Pseudo-Day'; colormap(jet(nanmax(pday)));
        case 'rday'
            scatter(dims(plt,1),dims(plt,2),[],rday(idx),'.'); str = 'Range of Pseudo-Days'; colormap(jet(nanmax(rday)));            
        case 'sday'
            scatter(dims(plt,1),dims(plt,2),[],sday(idx),'.'); str = 'Std of Pseudo-Days'; colormap(jet);
        case 'pdiff'
            scatter(dims(plt,1),dims(plt,2),[],pdiff(idx),'.'); str = '\DeltaPseudo-Day'; colormap(jet(nanmax(pdiff))); 
      
        % custom
        case 'custom'
            scatter(dims(plt,1),dims(plt,2),[],custom(idx),'.'); str = customstr; colormap(distinguishable_colors(length(unique(custom)),'w'))
        otherwise
            disp('invalid grouping!')
            return
    end    
    ti1 = title('Loading scatter plot...');
    
    xlabel([dimname ' 1'])
    c = colorbar;
    c.Label.String = str;
    ylabel([dimname ' 2'])
    
    %% select cluster
    if selectMode==1 || ~exist('ix','var')       
        ti1.String='Select points to visualize';
        if strcmpi(dimname,'tsne')
            dthresh = [5 5];
        else
            dthresh = diff(reshape(axis,2,2))/32;
        end
        h = imfreehand;
        pos = getPosition(h);
        in=inpolygon(dims(:,1),dims(:,2),pos(:,1),pos(:,2));
        in = in & plt;
        if ~any(in)
            in = abs(dims(:,1)-pos(1,1))<dthresh(1) & abs(dims(:,2)-pos(1,2))<dthresh(2);
            in = in & plt;
        end
        ix = kp(in);
        ti1.String='Loading points...';
    end
    if selectMode<2
        % three random trials
        if length(ix)>1
            ix2 = randsample(ix,min([3 length(ix)]));  ix2 = sort(ix2);
        else
            ix2 = ix;
        end
        
        plotted1=false; plotted2=false;

        selectMode = 2; % select trials
        value = zeros(1,3); r=1;
        loaded = false(size(ix2));
        mov = cell(length(ix2),1);
    end
    
    %% interactively select trials
    while selectMode==2
        
        %% heatmap plot        
        if ~plotted1 || ~plotted2
            h2=figure(2);
            WinOnTop(h2);
            set(h2,'pos',fpos(2,:))
            figure(2)
            plotMultiHeatmap(analysis,ix,ix2);            
            plotted1=true;
        end
        
        %% 2-D dom paw trajectories (three random trials)
        if ~plotted2
            figure(3)
            set(3,'pos',fpos(3,:))
            plot2DPaw(analysis, ix2, paw)
            
            [~,ix3] = ismember(ix2,kp);
            ti3 = title(sprintf('%d    ',round(pdist(trajPC(ix3,:))/mthresh*100)));
        end
        
        %% speed profiles (three random trials)
        if ~plotted2
            figure(4)
            set(4,'pos',fpos(4,:))
            plotSpeedProfile(analysis, ix2, ix)
        end
        plotted2 = true;
        
        %% update trials
        figure(2)
        if showVideo
            ti1.String = 'Select trials or press 0-3 to play a video';
        else
            ti1.String = 'Select trials or press any key';
        end
        
        saveFlag = false;
        k = waitforbuttonpress;
        if k==0 % mouse
            ti1.String = 'Loading trials...';
            value = round(get(gca,'CurrentPoint'));
            value(2,:) = [];
            if value(2)>=1 && value(2)<=length(ix)
                if select3
                    loaded = false(size(ix2));
                    ix2 = zeros(1,3);
                    while max(ix2)==0
                        try
                            n = round(length(ix)/30);
                            ix2 = ix(value(2)+[-randi(n) 0 randi(n)]);
                        end
                    end
                else
                    loaded(r) = false;
                    ix2(r) = ix(value(2));
                end
                plotted2 = false;
            end
        else
            value = double(get(gcf,'CurrentCharacter')); %48:51 --> 0 thru 3
            if value==110 % n = 2-nearest neighbors to blue trace
                ti1.String = 'Finding nearest neighbors...';
                knn = nnidx(ix2(1),:); knn(isnan(pdays(ix2(1),:)))=[];
                if length(knn)>=neighbors(2)
                    ix2(2:3) = knn(neighbors);
                    loaded(2:3) = false;
                elseif length(knn)>=neighbors(1)                    
                    ix2 = [ix2(1), knn(neighbors(1))];
                    loaded = [true false];
                else
                    ix2 = ix2(1);
                    loaded = true;                       
                end
                plotted2 = false;
            elseif value==115 % s = save
                ti1.String = 'Choose trial(s) to save';
                k = waitforbuttonpress;
                value = double(get(gcf,'CurrentCharacter')); %48:51 --> 0 thru 3
                if ~(value>=48 && value<=51)
                    value = 115;
                else
                    saveFlag = true;
                end
            end
        end
        
        %% videos
        
        if showVideo>0 && length(value)==1 && value>=48 && value<=51
            %% load videos
            ti1.String = 'Loading videos..';
            pause(0.01)
            
            for i = 1:length(ix2)
                if value==48 || value-48==i
                    if ~loaded(i)
                        % load video
                        mov{i} = playVideoTraj(analysis,ix2(i),1);
                        loaded(i) = true;
                    end
                    ti1.String = [ti1.String '..' num2str(i) '.'];
                    pause(0.01)
                end
            end
            
            for i = 1:length(ix2)
                if value==48 || value-48==i
                    if ~saveFlag
                        ti1.String = sprintf('Playing video %d',i);
                    else
                        ti1.String = sprintf('Saving video %d',i);
                    end
                    %% initial plot
                    figure(5)
                    clf
                    if ~saveFlag
                        set(5,'pos',fpos(5,:),'color','w')
                    else
                        set(5,'pos',fpos(6,:),'color','w')
                    end
                    playVideoTraj(analysis,ix2(i),mov{i},showVideo,saveFlag,i);
                end
            end
            value = zeros(1,3);
            %%
        elseif length(value)==1 && value==32 % space bar  
            ti1.String = 'Changing scatter plot...';
            g = mod(g,length(groups))+1;
            selectMode=2; % reset scatter but stay in select-trials mode
            break
        elseif length(value)==1 && value~=110 && value~= 115
            selectMode=1; % switch to select trials mode            
        end
        figure(1)
    end
    
end
try
WinOnTop(h1,false);
WinOnTop(h2,false);
end
out = {ix, ix2, [sessnum(ix) trialnum(ix)], [sessnum(ix2) trialnum(ix2)]};


return
% catch err
%    fprintf('Error on line %d\n',err.stack(end).line)
%    disp(err.message)
% end

