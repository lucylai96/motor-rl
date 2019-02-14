function results = twotap_agent(O,T,lesioned,plt,lrn,div)
% PURPOSE: simulate the agent navigating the 2-tap task, actor-critic learning for POMDPs
% AUTHOR: lucy lai
%
% USAGE: results = TD(x,O,T)
%
% INPUTS:
%   x - last observation (1 = nothing, 2 = tone)
%   a - last action (1 = nothing, 2 = tap)
%   O - [S' x A x O] observation distribution: O(i,j,m) = P(x'=m|s'=i,a=j)
%   probability of observing x in state s' after taking action a
%   T - [S x S' x A] transition distribution: T(i,j,k) = P(s'=j|s=i, a=k)
%   probability of transitioning to state s' from state s after taking
%   action a
%   lesioned = does the animal have MC or not
%
% OUTPUTS:
%   results - structure with the following fields:
%        .w - weight vector at each time point (prior to updating)
%        .b - belief state at each time point (prior to updating)
%        .rpe - TD error (after updating)
%

nTrials = div*10;
blur = 0.25;
C = 0.3; % cost to tap

%% lesioned?
if lesioned ==1
    T_lesion = T([1:12 25],[1:12 25],:);
    T_lesion(12,13,1) = 1;
    T_lesion(:,1,2) = 1;  % if you tap too early, state goes back to 1
    T_lesion(5,1,2) = 0.5; % except in these states
    T_lesion(6,1,2) = 0.5;
    T_lesion(7,1,2) = 0;
    T_lesion(8,1,2) = 0;
    T_lesion(7,13,2) = 1;
    T_lesion(8,13,2) = 1;
    T_lesion(13,13,2) = 0;
    
    for i = 1:size(T_lesion,1)
        G_blur_les(i,:) = normpdf(1:size(T_lesion,1),i,blur*i);
    end
    G_blur_les = [zeros(size(T_lesion,1),1) G_blur_les(:,1:12)];
    newT_lesion = T_lesion(:,:,1)+G_blur_les;
    newT_lesion=newT_lesion./sum(newT_lesion,2);
    T_lesion(:,:,1) = newT_lesion;
    
    O_lesion = O([1:12 25],:,:);
    
    O_lesion(13,2,1) = 0;
    O_lesion(13,2,2) = 1; %only way they see reward is when they have just tapped
end


for i = 1:25
    G_blur(i,:) = normpdf(1:25,i,blur*i);
end
G_blur = [zeros(25,1) G_blur(:,1:24)];

newT = T(:,:,1)+G_blur;
newT=newT./sum(newT,2);
T(:,:,1) = newT;

% figure;imagesc(T(:,:,1))
% set(gca,'YDir','normal')
% xlabel('s(t+1)')
% ylabel('s(t)')
% title('transition matrix: wait')


%% initialize states and weights
% initialization
S = size(T,1);      % number of states
b = ones(S,1)/S;    % belief state
if lesioned ==1
    b(13:25) = 0; %no belief in these states
end
%
% blur = 0.1; %proportionality constant for Gaussian blur
%
% for i = 1:S-1
%     G_blur(i,:) = normpdf(1:25,i,blur*i);
% end
% G_blur(25,:) = max(G_blur(:)).*ones(1,S);

%
% G_blur = G_blur./max(G_blur(:));
%b = b+G_blur(S0,:); %currently blurring belief by absolute time, may not be super plausible, maybe the reset happens when they tap
b = b/sum(b);

theta_policy  = zeros(S,2);     % policy weights
theta_policy(3,2) = 0.3;     % innate bias to tap at 300ms
theta_policy(4,2) = 0.3;     % innate bias to tap at 300ms

w_value  = zeros(S,1);          % value weights


%% initialize learning parameters
alpha_policy = 0.4;         % policy learning rate
alpha_value = 0.4;          % value learning rate
gamma = 0.98;
TE = 0.4; %temperature parameter in the policy (higher is more noisy)
last_S = S;


% figure; hold on;
% 
% subplot 121
% imagesc(T(:,:,1)); % for "null"
% title('a = null')
% axis square
% set(gca,'YDir','normal')
% 
% subplot 122
% imagesc(T(:,:,2)); % for "tone"
% axis square
% set(gca,'YDir','normal')
% title('a = tap')
% suptitle('transition matrices')




for t = 1:nTrials
    
    
    
    %% changing environment
    if lrn~=0
         T_new = twotap_changeworld(T,t,div,lrn);
         
        %if sum(T_new(:)-T(:)) ~=0
        % change(t) = 1;
       % end
         T = T_new;
    end
    
    %% take an action according to policy
    
    % if t ==1
    p_a = (exp((theta_policy'*b)./TE));%./(exp((theta_policy'*b)./T));
    %else
    %   p_action = exp(B*theta_policy'*results.b(t-1,:)');
    % end
    
    % if t==1
    %   p_a = [0 1]; %first thing agent will do is tap
    %else
    
    %end
    p_a=  p_a';
    
    %if t<10
    % p_a = [0.5 0.5];
    % a(t)=randi(2);
    
    %if t<50
    % if p_action(1) == p_action(2)
    %    p_a = [0.5 0.5];
    %   a(t)=1;
    % else
    
    p_a = p_a./sum(p_a);
    
    [~,a(t)] = max(p_a);
    
    %policy sampled stochastically
    if rand(1)<p_a(1)
        a(t) = 1;
    else
        a(t) = 2;
        
    end
    %end
    %else
    %    p_a = p_action./sum(p_action);
    %  [~,a(t)] = max(p_a);
    % end
    
    %% observe reward and new state
    [x(t), curr_S] = twotap_env(last_S,a(t), O, T);
    
    %if lesioned ==1
    %    if curr_S>12 && curr_S<25
    %        b0(curr_S) = 0;
    %        b0(curr_S-12) =1;
    %        b(curr_S) = 0;
    %        b(curr_S-12) =1;
    %    end
    %end
    
    %% belief state calculation
    
    if lesioned ==1
        b0= b; % old posterior, used later
        b = b([1:12,25]);
        b0_lesion= b; % old posterior, used later
        b = ((T_lesion(:,:,a(t))'*b0_lesion).*squeeze(O_lesion(:,a(t),x(t))));
        b = [b(1:12); zeros(12,1); b(13)];
        %[~,idx]=max(b');
        %b = b+G_blur(idx,:)'; %currently blurring belief by absolute time, may not be super plausible, maybe the reset happens when they tap
        %b = [b(1:12); zeros(12,1); b(13)];
        
    else
        
        b0 = b; % old posterior, used later
        %   [~,idx]=max(b0');
        %  b0 = b0+G_blur(idx,:)';
        % b0 = b0./sum(b0);
        b = ((T(:,:,a(t))'*b0).*squeeze(O(:,a(t),x(t))));
        
        %  [~,idx]=max(b');
        
        % She blurred the ISI distribution (see the first page of the online methods) and then used this blurred distribution to construct the transition matrix
        %  if t>1000
        % b = b+G_blur(idx,:)';
        %  end
        
        %b = conv(b,G_blur(idx,:)'); %currently blurring belief by absolute time, may not be super plausible, maybe the reset happens when they tap
        
    end
    
    
    
    %b=zeros(S,1);
    %b(curr_S) = 1;
    % b = b0'*(T(:,:,a(t)).*squeeze(O(:,a(t),x(t))));
    % b = b';
  
    b = b./sum(b); %normalize
    
    if any(isnan(b))
        keyboard
    end
    
    
    %% TD update
    w0 = w_value;
    r = double(x(t)==2)-C*(a(t)==2);        % reward (currently its binary)
    rpe = r + w_value'*(gamma*b-b0);  % TD error
    w_value = w_value + alpha_value*rpe.*b0;       % weight update
    
    % average reward replacing critic, (not state dependent)
    % critic not updating certain states
    
    
    
    %% policy update
    theta0 = theta_policy;      % reward
    
    %update thetas according to belief state-- theta(b(t))
    theta_policy(:,a(t)) = theta_policy(:,a(t)) + alpha_policy*gamma*rpe.*b0;   % policy weight update
    
    
    if any(isnan(theta_policy))
        keyboard
    end
    
    
    % gamma = gamma*gamma;
    %% store results
    results.t= t;
    results.w(t,:) = w0;
    results.b(t,:) = b0';
    results.theta(:,:,t) = theta0;
    results.rpe(t,:) = rpe;
    results.value(:,:,t) = w_value.*(b0); %estimated value
    results.l_state(t) = last_S;
    results.c_state(t) = curr_S;
    results.observe(t) = x(t);
    results.action(t) = a(t);
    results.s_action(t) = a(t)-1;
    results.p_action(t,:) = p_a;
    results.cost(t) = C;
    
    last_S = curr_S; %last state is now current state;
    
% if t>1
%         if a(t)==2 && a(t-1)==2
%             C = C+0.1;
%         %elseif a(t)==2
%         %    C = C+0.03;
%         elseif a(t)==2
%             C = C+0.2;
%         elseif C>0.03
%             C = C-0.03;
%             %theta_policy(25,2) = theta_policy(25,1)+.2;
%         end
% end

    
    if t>11
        if sum(results.s_action(t-11:t))>2 && C<0.4
            C = C+0.03;
        %elseif a(t)==2
        %    C = C+0.03;
        elseif C>0.1
            C = C-0.03;
            %theta_policy(25,2) = theta_policy(25,1)+.2;
        end
    end
    
    
    
end

%% plots
if plt ==1
    figure; hold on;
   % idx_change =find(change==1);
    
    % states over time
    subplot 511
    plot(results.l_state,'ro-');
    title('state')
    ylabel('state #')
   % line([idx_change' idx_change']',[repmat([0 26],10,1)]','Color','k')
    %xlabel('timesteps (a.u. ~100ms each)')
    
    % actions over time
    subplot 512
    plot(results.action,'go-');
    title('action chosen')
    ylabel(' action (1=wait, 2=tap)')
    %xlabel('timesteps (a.u. ~100ms each)')
    
    subplot 513
    imagesc(results.w');
    title('state value weights')
    ylabel('state #')
  %  line([idx_change' idx_change']',[repmat([0 26],10,1)]','Color','r')
    set(gca,'YDir','normal')
    %xlabel('timesteps (a.u. ~100ms each)')
    
    subplot 514
    imagesc(results.b');
    title('inferred belief state')
    ylabel('state #')
   % line([idx_change' idx_change']',[repmat([0 26],10,1)]','Color','r')
    set(gca,'YDir','normal')
    %xlabel('timesteps (a.u. ~100ms each)')
    
    subplot 515
    imagesc(squeeze(results.theta(:,2,:)));
    ylabel('state #')
    xlabel('timesteps (a.u. ~100ms each)')
    set(gca,'YDir','normal')
   % line([idx_change' idx_change']',[repmat([0 26],10,1)]','Color','r')
    title('policy weights')
    
    suptitle(strcat('total correct trials: ',num2str(sum(x>1))));
    
end

% figure
% for t = 1:length(x)
%     imagesc(results.b(t,:));
%     pause(.04);
% end


% value weights learned over time

% policy weights learned over time

end


