function [ks, metrics, stats, centroids_good] = loadKsMetrics(ksPath, channelMapFile, varargin)
%
%
% load basic kilosort metrics computed using neuropixel-utils

p = inputParser();
p.addRequired('ksPath',@isstringlike);
p.addRequired('mapfile', @isstringlike);
p.addParameter('saveResults',true, @islogical);


p.KeepUnmatched = false;
p.parse(ksPath, channelMapFile, varargin{:})

% /input parsing

channelMap = Neuropixel.ChannelMap(channelMapFile);
ks = Neuropixel.KilosortDataset(ksPath, 'channelMap', channelMap, deduplicate_spikes=false, deduplicate_cutoff_spikes=false);
ks.load(loadFeatures=false, loadBatchwise=true);

metrics = ks.computeMetrics();
stats = ks.computeBasicStats();
ks.printBasicStats();

% pull out centroid positions
metrics = ks.computeMetrics();

mask_good = ismember(metrics.cluster_ids, ks.clusters_good);
centroids_good = metrics.cluster_centroid(mask_good,:);


if p.Results.saveResults
    [~,recordingName,~] = fileparts(ksPath)
    filename = fullfile(ksPath,[recordingName '_ksMetrics.mat']);
    save(filename,'metrics','stats','ks','centroids_good','mask_good')

    % export centroid positions as .tsv
    header = {'cluster_id','centroid_x','centroid_y'};
    data = [double(ks.clusters_good), centroids_good];
    data = mat2cell(data, ones(size(data,1),1), ones(size(data,2),1));
    C = [header; data];

    filename = fullfile(ksPath,[recordingName '_cluster_centroids.tsv']);
    writecell(C,filename, 'filetype','text', 'delimiter','\t')

end