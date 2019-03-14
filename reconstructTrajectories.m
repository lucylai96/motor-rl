function analysis = reconstructTrajectories(analysis)
struct_unzip(analysis);

fprintf('reconstructing... '); tic
% un-rotate 
scores = trajPC*coeff2';

% extract each signal from PCs
data_rc = nan(data_dims);

ncum = [0 cumsum(ncomps)];
for s = 1:length(signals)
    
    i = []; j = 1:2;    
    switch signals{s}
        case 'paw'
            if ~any(strcmp(signals,'arm'))
                i = paw;
            end
        case 'paw2'
            i=paw2;
        case 'arm'
            i = [paw, elbow];
        case 'arm2'
            i = [paw2 elbow2];
        case 'head'
            i = head;
        case 'lever'
            i = anFlag; j=1;
        otherwise
    end
    if isempty(i) || i(1)==0
        fprintf('no %s!\n', signals{s})
    else 
        temp = scores(:,ncum(s)+1:ncum(s+1))*coeff{s} + pc_mu{s};
        data_rc(kp,:,j,i,:) = reshape(temp, size(data_rc(kp,:,j,i,:)));
    end
end

% reconstruct position and velocity
if any(strcmp(features,'position'))
    i = find(strcmp(features,'position'));
    if tapFlag
        mu = tap_mu; sd = tap_std;
    else
        mu = traj_mu; sd = traj_std;
    end
    temp = data_rc(:,:,:,:,i);
    temp = temp.*permute(sd,[2 3 4 1]);
    pos_rc = temp+mu;
end
if any(strcmp(features,'velocity'))
    i = find(strcmp(features,'velocity'));
    if tapFlag
        mu = tap_mu_vel; sd = tap_std_vel;
    else
        mu = traj_mu_vel; sd = traj_std_vel;
    end
    temp = data_rc(:,:,:,:,i);
    temp = temp.*permute(sd,[2 3 4 1]);
    vel_rc = temp+mu;
else
    if tapFlag
        t = diff(time_tap);
    else
        t = repmat(time_warp,length(warpFactor),1);
        j = find(time_warp>=0 & time_warp<=700);
        t(:,j) = warpFactor*time_warp(j);
        j = max(j)+1;
        t(:,j:end) = t(:,j:end)-(t(:,j)-t(:,j-1))+median(diff(time_warp));
        t = diff(t,[],2);
    end
    vel_rc = diff(pos_rc,[],2) ./ t;
    vel_rc = vel_rc(:,[1 1:end],:,:);
end
if ~exist('pos_rc','var')
    if tapFlag
        t = median(diff(time_tap));
        mu = tap0;
    else
        t = repmat(time_warp,length(warpFactor),1);
        j = find(time_warp>=0 & time_warp<=700);
        t(:,j) = warpFactor*time_warp(j);
        j = max(j)+1;
        t(:,j:end) = t(:,j:end)-(t(:,j)-t(:,j-1))+median(diff(time_warp));
        t = diff(t,[],2);
        mu = traj0;
    end
    pos_rc = nan(size(vel_rc)); pos_rc(:,1,:,:) = mu;
    pos_rc(:,2:end,:,:) = mu+cumsum(vel_rc(:,2:end,:,:).*t,2,'omitnan');
end

if tapFlag
    analysis.tap_filt_rc = pos_rc;
    analysis.tap_vel_rc = vel_rc;
else
    analysis.traj_warp_rc = pos_rc;
    analysis.traj_warp_vel_rc = vel_rc;
end
toc
return
