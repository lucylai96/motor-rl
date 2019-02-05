function [results_intact,results_lesioned]=twotap_world(plt)
% purpose: simulate the world in the 2-tap task
% macro-states = {1:start, 2:IPI, 3:ITI, 4:reward}
% micro-states = {1-13 is IPI, 14-25 is ITI, 26 is rew, 27 is start}

% action = {1:wait, 2:tap}
% observations = {1:null, 2:tone}


%% initialize
% can run the exp on two timescales: one at the level of second matter, one
% where at the level of trials
nSt= 26;%number of states
nAc = 2; %number of actions
nOb= 2;%number of observations


%% fill out the observation matrix O

%O - [S' x O x A] observation distribution: O(i,j,m) = P(x'=m|s'=i,a=j)
%   probability of observing x in state s' after taking action a

O=zeros(nSt,nAc,nOb); % 26 states, 3 actions, 2 observations


O(26,2,2) = 1; % only state-action pair where you observe tone is tap (#2) )
%O(26,2,2) = 1; % only state-action pair where you observe tone is tap (#2)

O(:,1,1) = 1; % in states 1-25,27, you will observe nothing if you wait (#1)
O(:,2,1) = 1; % in all states, if you tap (#2) you will observe nothing
O(26,2,1) = 0;
%O(26,2,2) = 0;

figure; hold on;
subplot 121
imagesc(O(:,:,1)); % for "null"
title('x = null')
subplot 122
imagesc(O(:,:,2)); % for "tone"
title('x = tone')
suptitle('observation matrices')

%% fill out the transition matrix T
%T(i,j,k) is the probability of transitioning from sub-state i-->j after
%taking action k
T= zeros(nSt,nSt, nAc);

% page 1: wait
T(1,2,1)=1;
T(2,3,1)=1;
T(3,4,1)=1;
T(4,5,1)=1;
T(5,6,1)=1;
T(6,7,1)=1;
T(7,8,1)=1;
T(8,9,1)=1;
T(9,10,1)=1;
T(10,11,1)=1;
T(11,12,1)=1;
T(12,13,1)=1;
T(13,26,1)=1;

T(14,15,1)=1;
T(15,16,1)=1;
T(16,17,1)=1;
T(17,18,1)=1;
T(18,19,1)=1;
T(19,20,1)=1;
T(20,21,1)=1;
T(21,22,1)=1;
T(22,23,1)=1;
T(23,24,1)=1;
T(24,25,1)=1;
T(25,26,1)=1;
%T(26,27,1)=1; %if you wait you go to 27
T(26,26,1)=1;

% page 2: tap
T(1,14,2)=1;
T(2,14,2)=1;
T(3,14,2)=1;
T(4,14,2)=1;
T(5,14,2)=0.5;
T(6,14,2)=0.5;

T(5,26,2)=0.5;
T(6,26,2)=0.5;
T(7,26,2)=1;
T(8,26,2)=1;

%T(8,14,2)=1;
T(9,14,2)=1;
T(10,14,2)=1;
T(11,14,2)=1;
T(12,14,2)=1;
T(13,14,2)=1;

T(14,14,2)=1;
T(15,14,2)=1;
T(16,14,2)=1;
T(17,14,2)=1;
T(18,14,2)=1;
T(19,14,2)=1;
T(20,14,2)=1;
T(21,14,2)=1;
T(22,14,2)=1;
T(23,14,2)=1;
T(24,14,2)=1;
T(25,14,2)=1;
%T(26,27,2)=1; %if you tap you go to 27 (redundant)
T(26,1,2)=1;

%figure; hold on
%subplot 121
%imagesc(T(:,:,1)); % for "nothing"
%title('A = nothing')
%subplot 122
%imagesc(T(:,:,2)); % for "tap"
%title('A = tap')
%suptitle('transition matrices')




%% "intact" animal

results_intact= twotap_agent(O,T,0,plt);


%% "lesioned" animal

results_lesioned = twotap_agent(O,T,1,plt);
% the lesioned animal only knows the IPI and not the ITI
% wait, so the transition matrix should stay the same, but

end