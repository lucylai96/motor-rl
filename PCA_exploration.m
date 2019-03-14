function PCA_exploration
close all
% purpose: this function examines how the transition between different taps
% as defined in PCA space is

% a goal is to be able to visualize how taps are explored in a space for% a goal is to be able to visualize how taps are explored in a space for
% lesioned vs nonlesioned animals


%% adding dependencies
addpath('/Users/lucy/Google Drive File Stream/Team Drives/MC Learning Project/Matlab/output/')
addpath('/Users/lucy/Google Drive File Stream/Team Drives/MC Learning Project/Matlab/gerald new code/')

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


% 16889 trials
tap_exploration('GP1829_001_061_taps_v_aah.mat',1) %naive partial lesion but seemed to learn IPI kind of (after 5000 trials), and ITI kind of!!

% 16869 trials
tap_exploration('GP1830_001_084_taps_v_aah.mat',2)% this guy was switched to auditory feedback but before then did not learn ANY structure

% 20321 trials
tap_exploration('GP1831_001_045_taps_v_aah.mat',3)% this guy systematically undershot IPI, no stereotyped ITI but did learn some low CV for IPI

% 17084 trials
tap_exploration('GP1832_001_070_taps_v_aah.mat',4) %learned around 13k trials

% 16751 trials
tap_exploration('GP1838_001_099_taps_v_aah.mat',5)% learned around 5k trials

% 17800 trials
tap_exploration('GP1840_001_059_taps_v_aah.mat',6)% learned around 5k trials
end

function tap_exploration(rat,s)



%% here add different resolutions in ways of parsing the data
%% plot distance against time
analysis = load(rat);


%% separate by first and second tap in a trial

%% separate by
end