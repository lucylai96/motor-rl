function results = twotap_agent(S0,O,T,lesioned)
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

nTrials = 3000;

%% initialize states and weights
% initialization
S = size(T,1);      % number of states
b = ones(S,1)/S;    % belief state
%b(26) = 1;
b = b/sum(b);
theta_policy  = zeros(S,2);     % policy weights
%theta_policy(:,1) = 0.2;     % policy weights
theta_policy(end,2) = 0; % always tap at the beginning
theta_policy([7,26],2) = 0.2;
%theta_policy(14,1) = 1;
%theta_policy([1:6,8:13,14:25],1) = 0.7;
w_value  = zeros(S,1);          % value weights
C = 0.2;%cost to tap

%% initialize learning parameters
alpha_policy = 0.2;         % policy learning rate
alpha_value = 0.2;          % value learning rate
gamma = 0.98;
TE = .3; %temperature parameter in the policy (higher is more noisy)
last_S = S;

for t = 1:nTrials
    
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
    
    if lesioned ==1
        if curr_S>13 && curr_S<26
            b(curr_S) = 0;
            b(curr_S-13) =1;
        end
    end
    
    %% belief state calculation
    b0 = b; % old posterior, used later
    b = ((T(:,:,a(t))'*b0).*squeeze(O(:,a(t),x(t))));
    %b=zeros(S,1);
    %b(curr_S) = 1;
    % b = b0'*(T(:,:,a(t)).*squeeze(O(:,a(t),x(t))));
    % b = b';
    
    b = b./sum(b); %normalize
    % b(14:25) = 0;
    
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
    if any(isnan(b))
        keyboard
    end
    
    % gamma = gamma*gamma;
    %% store results
    results.w(t,:) = w0;
    results.b(t,:) = b0';
    results.theta(:,:,t) = theta0;
    results.rpe(t,:) = rpe;
    results.value(:,:,t) = w_value.*(b0); %estimated value
    results.state(t) = last_S;
    results.observe(t) = x(t);
    results.action(t) = a(t);
    results.p_action(t,:) = p_a;
    results.cost(t) = C;
    
    last_S = curr_S; %last state is now current state;
    
    if C <0.3
        C = C+0.03;
    else
        C = C-0.05;
        %theta_policy(26,2) = theta_policy(26,1)+.2;
    end
    
    
    
end

%% plots

figure; hold on;

% states over time
subplot 411
plot(results.state,'ro-');
title('state')
ylabel('state #')
%xlabel('timesteps (a.u. ~100ms each)')

subplot 412
imagesc(results.w');
title('state value weights')
ylabel('state #')
%xlabel('timesteps (a.u. ~100ms each)')

subplot 413
imagesc(results.b');
title('inferred belief state')
ylabel('state #')
%xlabel('timesteps (a.u. ~100ms each)')

subplot 414
plot(results.p_action(:,2)');
ylabel('p(tap)')
xlabel('timesteps (a.u. ~100ms each)')
title('probability of tapping')

suptitle(strcat('total correct trials: ',num2str(sum(x>1))));
%figure
%plot(results.cost,'bo-');
% value weights learned over time
for t = 1:length(x)
    results.w(t,:)
end
% policy weights learned over time

figure(1);hold on;
%plots
temp = exp(results.theta./TE);
temp = temp./sum(temp,2);
for i = 1:length(x)
    
    K(i,:) = temp(:,1,i)';
    plot(K(i,:))
    pause(0.01);
end
imagesc(K);
end

% function [x,next_S] = twotap_env(last_S,a, O, T)
% % takes in action and previous state, outputs observation x
%
% [~, next_S] = max(T(last_S,:,a));
% [~, x] = max(O(next_S,a,:));
%
% end
