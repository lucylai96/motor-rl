function [x,next_S] = twotap_env(last_S,a, O, T)
% takes in action and previous state, outputs observation x

[~, next_S] = max(T(last_S,:,a));
[~, x] = max(O(next_S,a,:));

end