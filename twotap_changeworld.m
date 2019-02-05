function T = twotap_changeworld(T, tNum,div)
% inputs:
%      T_tap = second "page" of transition matrix, what happens when you tap.
%      Basically assigns which intervals get reward.
%      tNum = trial number

%div = 1000;


% if mod(tNum,div) ==1%
%     figure(100); hold on;
%     subplot 121
%     imagesc(T(:,:,1)); % for "nothing"
%     set(gca,'YDir','normal')
%     title('A = nothing')
%     
%     subplot 122
%     imagesc(T(:,:,2)); % for "tap"
%     title('A = tap')
%     set(gca,'YDir','normal')
%     suptitle('transition matrices')
%    
% end


if tNum>div && tNum<div*2
    
    T(3,14,2)=1;
    T(3,26,2)=0;
    
    T(5,26,2)=1;
    T(5,14,2)=0;
    
elseif tNum>div*2 && tNum<div*3
    %T(7,26,2)=1;
    %T(8,26,2)=1;
    
    T(4,14,2)=1;
    T(4,26,2)=0;
    
    T(6,26,2)=1;
    T(6,14,2)=0;
    
elseif tNum>div*3 && tNum<div*4
    
    T(5,14,2)=1;
    T(5,26,2)=0;
    
    T(7,26,2)=1;
    T(7,14,2)=0;
    
elseif tNum>div*4 && tNum<div*5
    
    T(6,14,2)=1;
    T(6,26,2)=0;
    
    T(8,26,2)=1;
    T(8,14,2)=0;
    
    %% seven-five-seven
elseif tNum>div*5 && tNum<div*6
    

    
elseif tNum>div*6 && tNum<div*8
    
    T(7,14,2)=1;
    T(7,26,2)=0;
    T(8,14,2)=1;
    T(8,26,2)=0;
    
    T(5,26,2)=1;
    T(5,14,2)=0;
    
elseif tNum>div*8 && tNum<div*10
      T(5,14,2)=1;
         T(5,26,2)=0;
         
            T(7,26,2)=1;
    T(7,14,2)=0;
    
end




end