function [spikeIdxMat, clusterId, clusterLabels] = ksResults2spikeMat(ksResultsPath)
% Loads spike times from a kilosort results folder and packs them into a [nSample x nUnit] sparse array 
% 
% *** Note: spike times are in units of samples, not time ***
%
% EMT 2021-03-22
%
%
% #TODO: add support for loading different cluster groups (good, mua,unsorted, noise)


tSpk = readNPY(fullfile(ksResultsPath,'spike_times.npy'));
clusterId = readNPY(fullfile(ksResultsPath,'spike_clusters.npy'));
label = tdfread(fullfile(ksResultsPath,'cluster_group.tsv'),'\t');

label.group = mat2cell(label.group, ones(size(label.group,1),1), size(label.group,2));
label.group = cellfun(@(x) strtrim(x), label.group, 'uni',false);

goodUnit = cellfun(@(x) strcmp(x,'good'), label.group);
multiUnit = cellfun(@(x) strcmp(x,'mua'), label.group);
noise = cellfun(@(x) strcmp(x,'noise'), label.group);
% nUnit = nnz(goodUnit | multiUnit );
nSample = max(tSpk);    % find the index of the last spike in the dataset as a proxy for the length of the recording 

% group spikes by unit
spkIdx = cellfun(@(x) tSpk(clusterId==x),...
    num2cell(label.cluster_id(goodUnit)), 'uni',false);
warning('ignoring all multiunits - todo: add functionality to keep these')

spkIdx = cellfun(@unique, spkIdx, 'uni', false);
% spkIdx = cellfun(@(x) double(x), spkIdx, 'uni',false);

% pack [nSample x nUnit] sparse array of spike times 
spikeIdxMat = cell2mat(cellfun(@(idx) sparse(double(idx),1,true,double(nSample),1), spkIdx', 'uni',false));

% for cluster labels
clusterMask = cellfun(@(x) strcmp(x,'good') | strcmp(x,'mua'), label.group);

clusterLabels = label;
clusterLabels.cluster_id = clusterLabels.cluster_id(clusterMask);
clusterLabels.group = clusterLabels.group(clusterMask);
