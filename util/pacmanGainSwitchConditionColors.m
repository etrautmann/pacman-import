function [colors] = pacmanGainSwitchConditionColors(Conditions)
%
%
% EMT 2021-03-25


gain = Conditions.gain;

colors = zeros(size(Conditions,1), 3);

colors(gain==1,:) = repmat([1 0 0],nnz(gain== 1),1);
colors(gain== -1,:) = repmat([0 0 1],nnz(gain== -1),1);
