function [TS] = exportTaskTable(Task, fileName, varargin)
%
% Work in progress
%
% EMT / HSC 2021-04-01

% Conditions
conds = table2struct(Task.Conditions);

% Target force
targetForces = makerow(Task.targetForce);

% Beh forces
forces = {Task.Force.data};

% Beh force labels
forceLabels = {Task.Force.variableLabels};
forceLabels = forceLabels{1};

% Alignment indices
alignIdxs = {Task.Force.alignIndex};

% Spikes
spikes = {Task.MU.spikes};

% Single trial
rates = {Task.MU.rate};

% PSTH (trial average)
psths = cellfun(@(x) squeeze(x(:,:,1)), {Task.MU.psth}, 'UniformOutput', false);    % keep mean, strip variance


% mask out specific conditions:
condMask = cellfun(@(x) size(x,3) > 1, forces);
    
conds = conds(condMask);
targetForces = targetForces(condMask);
forces = forces(condMask);
alignIdxs = alignIdxs(condMask);
spikes = spikes(condMask);
rates = rates(condMask);
psths = psths(condMask);

n_conds = length(conds);

sprintf("saving data to: %s",fileName)
save(fileName, 'n_conds','conds','targetForces','forces','forceLabels','alignIdxs','spikes','rates','psths', '-v7.3')

