% need trajTSNE as input

function [bins]=watershed_sandbox(trajTSNE,s)
%% generate density /watershed
% close all
maxVal = max(abs(trajTSNE(:)));
maxVal = round(maxVal * 1.1);
sigma = maxVal / s;
numPoints = 1001;
rangeVals = [-maxVal maxVal];

%figure
[xx,density] = findPointDensity(trajTSNE,sigma,numPoints,rangeVals);
%imagesc(xx,xx,density)
%set(gca,'ydir','normal')
%hold on
L = watershed(-density,8);

[ii,jj] = find(L==0);
%plot(xx(jj),xx(ii),'k.')

% figure
% imagesc(xx,xx,L)
% set(gca,'ydir','normal')
% cmap = [1 1 1; distinguishable_colors(max(L(:)))];
% colormap(cmap)

%% get watershed bin
vals = round((trajTSNE + max(xx))*length(xx)/(2*max(xx)));

N = size(trajTSNE,1);
watershedMode = zeros(N,1);
for i = 1:N    
    watershedMode(i) = L(vals(i,2),vals(i,1));
    if watershedMode(i) == 0        
        temp = L(vals(i,2)+(-1:1),vals(i,1)+(-1:1));
        watershedMode(i) = mode(temp(temp>0));
    end
end
n=hist(watershedMode,1:max(watershedMode));

goodmodes = find(n>10);

watershedMode(~ismember(watershedMode,goodmodes)) = 0;
[~,~,watershedMode] = unique(watershedMode);

n = max(watershedMode);
cmap = parula(n);

bins = watershedMode;


 scatter(trajTSNE(:,1),trajTSNE(:,2),[],watershedMode,'.')
 colormap(cmap)
 c = colorbar;
 c.Label.String = 'Mode';
 ticks = linspace(1,n,n*2+1); ticks = ticks(2:2:end);
 c.Ticks = ticks;
 c.TickLabels = 1:n;
 ylabel('tSNE 2')
 xlabel('tSNE 1')
end