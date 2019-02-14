function T = twotap_changeworld(T, tNum,div,lrn)
% inputs:
%      T_tap = second "page" of transition matrix, what happens when you tap.
%      Basically assigns which intervals get reward.
%      tNum = trial number

%div = 1000;


if mod(tNum,div) ==1%
    figure(100); hold on;
    subplot 121
    imagesc(T(:,:,1)); % for "nothing"
    set(gca,'YDir','normal')
    title('A = nothing')

    subplot 122
    imagesc(T(:,:,2)); % for "tap"
    title('A = tap')
    set(gca,'YDir','normal')
    suptitle('transition matrices')
   % pause(.1)

end

if lrn ==1

if tNum>div && tNum<div*2
    
    T(3,13,2)=1;
    T(3,25,2)=0;
    
    %  T(5,25,2)=1;
    % T(5,13,2)=0;
    
elseif tNum>div*2 && tNum<div*3
    %T(7,25,2)=1;
    %T(8,25,2)=1;
    
    T(4,13,2)=1;
    T(4,25,2)=0;
    
    %T(6,25,2)=1;
    %T(6,13,2)=0;
    
elseif tNum>div*3 && tNum<div*6
    
    T(5,13,2)=1;
    T(5,25,2)=0;
    
    %T(7,25,2)=1;
    %T(7,13,2)=0;
    
elseif tNum>div*6 && tNum<div*10
    
    T(6,13,2)=1;
    T(6,25,2)=0;
    
    T(8,25,2)=1;
    T(8,13,2)=0;
    
end
    %% seven-five-seven
if lrn==2
    if tNum>0  && tNum<div*3
    T(6,13,2)=0;
    T(6,25,2)=1;
    
    T(7,13,2)=0;
    T(7,25,2)=1;
    
    T(8,25,2)=1;
    T(8,13,2)=0;
    
    
    
elseif tNum>div*3 && tNum<div*7
    
    T(7,13,2)=1;
    T(7,25,2)=0;
    T(8,13,2)=1;
    T(8,25,2)=0;
    
    T(5,25,2)=1;
    T(5,13,2)=0;
    T(4,13,2)=0;
    T(4,25,2)=1;
    
elseif tNum>div*7 && tNum<div*10
    T(5,13,2)=1;
    T(5,25,2)=0;
    T(4,13,2)=1;
    T(4,25,2)=0;
    
    
    T(7,25,2)=1;
    T(7,13,2)=0;
    T(8,25,2)=1;
    T(8,13,2)=0;
    
    end
end




end