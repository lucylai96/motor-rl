function [paramDensity, paramMean, paramCV] = plotDensityLearning(paramAll)
% simplified version of convParam
% plot a smoothed density of IPI (or ITI)

% input: a column vector of IPIs
numTrials = length(paramAll);

%% parameters
param = 'IPI'; % choose IPI or ITI

% start/stop are row vectors that can be used to prevent smoothing over
% change boundaries (e.g. boundary from 700 target to 500)
start = 1;
stop = numTrials;

rewardBounds = [500 700 900]; % lower, target, upper

plotFlag = true;

mwin = 300; % moving window size (trials)
if strcmp(param,'IPI')
    paramBins = 0:1200; % discretization of IPI
    upperScale = 7E-3; % max intensity in plot
    trialBinSize = 50; % blur trials
    intertapBinSize = 50; % blur IPI bins    
else
    paramBins = 0:5000;
    upperScale = 2E-3;
    trialBinSize = 100;
    intertapBinSize = 75;
end
%% compute density
bins = discretize(paramAll,paramBins);
inc = find(~isnan(bins));

binaryMatrix = full(sparse(bins(inc),inc,1,length(paramBins),numTrials)');
convMatrix = ones(trialBinSize, intertapBinSize);
paramDensity = nan(numTrials,length(paramBins));
for i = 1:length(start)
    temp = conv2(binaryMatrix(start(i):stop(i),:), convMatrix, 'same');
    paramDensity(start(i):stop(i),:) = temp;
end
paramDensity = bsxfun(@rdivide,paramDensity,sum(paramDensity,2));

%% compute mean/std/cv
paramMean = [];
paramMed = []; % use median for ITI
paramStd = [];

for i = 1:length(start)
    temp = paramAll(start(i):stop(i));
    paramMean = cat(2,paramMean,movmean(temp,mwin,'omitnan'));
    paramMed = cat(2,paramMed,movmedian(temp,mwin,'omitnan'));
    paramStd = cat(2,paramStd,movstd(temp,mwin,'omitnan'));
end
paramCV = paramStd./paramMean;

%% generate plot
if plotFlag
    paramDensity(paramDensity<0)=0;
    paramDensity(paramDensity>upperScale)=upperScale;
    paramDensity = paramDensity/upperScale;
    X = gray2ind(paramDensity');
    
    if strcmp(param,'IPI')
        yyaxis(gca,'left')
        set(gca,'ycolor','k')
    end
    
    % plot density
    image(X)
    colormap(gray)
    set(gca,'ydir','normal')
    hold on
    
    if strcmp(param,'IPI')
        % plot reward bounds
        plot(xlim,rewardBounds([1 1]),'-m')
        plot(xlim,rewardBounds([3 3]),'-m')
        plot(xlim,rewardBounds([2 2]),'-m','linewidth',1.5)
        ylim([0 1200])
        set(gca,'ytick',[100:200:1100])

        % plot mean
        plot(paramMean,'g','linewidth',2)
        ylabel('IPI')

        % plot CV on different axis
        yyaxis(gca,'right')
        set(gca,'ycolor',lines(1))
        plot(xlim,[0.25 0.25],':','color',lines(1))
        plot(paramCV,'-','color',lines(1),'linewidth',2)
        ylim([0 1])
        ylabel('CV')
    else
        plot(1200*ones(size(paramMean)),'m','linewidth',1.5)
        % plot median
        plot(paramMed,'g','linewidth',2)
        ylabel('ITI')
        ylim([0 2400])
    end
end