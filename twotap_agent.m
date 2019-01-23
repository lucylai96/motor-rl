function results = twotap_agent(S0,O,T)
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
%
% OUTPUTS:
%   results - structure with the following fields:
%        .w - weight vector at each time point (prior to updating)
%        .b - belief state at each time point (prior to updating)
%        .rpe - TD error (after updating)
%

nTrials = 10000;

%% initialize states and weights
% initialization
S = size(T,1);      % number of states
b = ones(S,2)/S;    % belief state
theta_policy  = zeros(S,1);     % policy weights
w_value  = zeros(S,1);          % value weights

%% initialize learning parameters
alpha_policy = 0.1;         % policy learning rate
alpha_value = 0.1;          % value learning rate
gamma = 0.98;
B = 1; %beta??
last_S = S0;

for t = 1:nTrials
    
    %% take an action according to policy
    if t ==1
    p_action = exp(B*theta_policy'*b);
    else
        p_action = exp(B*theta_policy'*results.b(:,:,t-1));
    end
    
    if p_action(1) == p_action(2)
        p_a = 0.5;
        if rand() >0.4;
        a(t) = 2;
        else
         a(t) = 1;
        end
    else
        p_a = p_action./sum(p_action);
        [~,a(t)] = max(p_a);
    end
    
    %% observe reward and new state
    [x(t), curr_S] = twotap_env(last_S,a(t), O, T);
    
    
    %% belief state calculation
     if t ==1
    b0 = b; % old posterior, used later
     else
         b0 = results.b(:,:,t-1); % old posterior, used later
     end
    b = ((T(:,:,a(t))*b0(:,a(t))).*squeeze(O(:,a(t),x(t))));
    
    %(T*b2)'*squeeze(O(:,:,x(t)))
    b = b./sum(b); %normalize
    
    
    
    %% TD update
    w0 = w_value;
    r = double(x(t)==2);        % reward
    rpe = r + w_value'*(gamma*b-b0(:,a(t)));  % TD error
    w_value = w_value + alpha_value*rpe.*b0(:,a(t));       % weight update
    
    %% policy update
    theta0 = theta_policy;      % reward
    theta_policy = theta_policy + alpha_policy*B*rpe*log(p_a);         % weight update
    
    %% store results
    results.w(t,:) = w0;
    results.b(:,:,t) = b0;
    results.theta(t,:) = theta0;
    results.rpe(t,:) = rpe;
    results.value(:,:,t) = w_value.*(b0); %estimated value
    results.state(t) = last_S;
    results.observe(t) = x(t);
    
    last_S = curr_S; %last state is now current state;
end

end
