function [Task, Tsync, ks, metrics, stats, centroids_good ] = import_pacman_data(subject,date, varargin)
% import raw neural and behavioral data for pacman task / Rig J
% EMT 2023-02-13
%
% 2023-02-13 Note: this script is intended to be a function that replaces the
% functionality of the script run_load_neural_and_behavior.m, so we can
% call this from python or more simply automate processing of multiple
% datasets.
%
% dependencies:
%   spikes:             https://github.com/cortex-lab/spikes
%   sortingQuality:     https://github.com/cortex-lab/sortingQuality
%   neuropixel-utils:   https://github.com/djoshea/neuropixel-utils
%   matlab-utils:       https://github.com/djoshea/matlab-utils
%   npy-matlab:         https://github.com/kwikteam/npy-matlab
%

p = inputParser();
p.addRequired('subject', @isstringlike);
p.addRequired('date', @isstringlike);
p.addParameter('gNum',0, @isnumeric);
p.addParameter('tNum',0, @isnumeric);
p.addParameter('imecNums',0, @isnumeric);
p.addParameter('padDuration', [.5, 1], @isnumeric);
p.addParameter('saveTags',0:20, @isnumeric);
p.addParameter('syncChan',3, @isnumeric);
p.addParameter('syncBit',8, @isnumeric);
p.addParameter('channelMapFile','neuropixNHPv1_kilosortChanMap_v1.mat', @isstringlike);

p.addParameter('loadNeural',true, @islogical);
p.addParameter('loadMetrics',true, @islogical);
p.addParameter('alignTrials', false, @islogical);
p.addParameter('saveTaskTable',true, @islogical);
p.addParameter('runConnectivityAnalysis',true,@islogical)
p.addParameter('saveResults', true, @islogical);
% p.addParameter('use_cached', true, @islogical);     % if analysis has already run with saved results, load these instead of recomputing

p.KeepUnmatched = false;    
p.parse(subject, date, varargin{:});

subject = p.Results.subject;
date = p.Results.date;
gNum = p.Results.gNum;
tNum = p.Results.tNum;
imecNums = p.Results.imecNums;

% /input parsing


 % check that environment variables are set
assert(~strcmp(getenv('DATA_ROOT'),''),'environment variable "DATA_ROOT" must be set')
assert(~strcmp(getenv('FIG_ROOT'),''),'environment variable "DATA_ROOT" must be set')

figPath = fullfile(getenv('FIG_ROOT'),'pacman-gain-switch',date);
mkdir(figPath)

% 0.3) setup paths
paths = pacmanPaths(getenv('DATA_ROOT'), subject, date, gNum, tNum, imecNums);  % TODO: add support for multiple probes in pacmanPaths

% manually override taskTableOutputPath (if you want to specify where to store processed task data)
% paths.taskTableOutputPath = fullfile(getenv('DATA_ROOT'), subject, 'processed', date, 'mergedTaskData', [paths.prefix '_taskdata.mat']);
% paths.taskTableOutputPath = paths.recordingRoot; 
% makeContainingFolder(paths.taskTableOutputPath)


%% 1) load behavior data

tic

T = loadsession(fullfile(paths.sgDataPath, paths.prefixBehavior));
fprintf('Session loaded: %d trials \n',size(T,1))


% attempt to load EMG data
% nsx = openNSx(paths.brDataPath,'uv');

timing = [];
timing.behaviorLoad = toc


%% 2) load kilosort output


if (p.Results.loadNeural == 1)
    [spikeIdxMat, clusterID, clusterLabels] = ksResults2spikeMat(paths.ksResultsPath);
    fprintf('Spike Times loaded %.1f \n', toc)
    spikeIdxMat_orig = spikeIdxMat;
    
    timing.spikeMatCreated = toc
    
    % 3) convert spiketimes from imec headstage sample indices into NIDAQ indices:
    niMeta = readSpikeGLXmeta(paths.nidaqMetaPath);
    FsNi = niMeta.sRateHz;
    
    apMeta = readSpikeGLXmeta(paths.npixApMetaPath);
    FsImec = apMeta.sRateHz;
    
    spikeIdxMatNi = convertSpikeTimeIndices(spikeIdxMat, FsImec, FsNi);
    size(spikeIdxMatNi)
    
    timing.spikeTimesConverted = toc
    
    % 4) sync spike times with behavior
    [Tsync, Times] = syncSpeedgoatNidaq(paths.nidaqPath, T, 'spikes',spikeIdxMatNi, 'SGsyncChan',p.Results.syncChan, 'SGsyncBit', p.Results.syncBit);
    fprintf('Spike Times synced: %.1f sec\n',toc)
    
    timing.nidaqSynced = toc
else
    Tsync = T; 
end

%% load neural metrics
if p.Results.loadMetrics
    [ks, metrics, stats, centroids_good] = loadKsMetrics(paths.ksResultsPath, p.Results.channelMapFile);
else
    ks = [];
    metrics = [];
    stats = [];
    centroids_good = [];
end


%% 5) Calculate condition info 


Task = paccond_gain_switch(Tsync,'neuropixels','saveTags',p.Results.saveTags,'padDur',p.Results.padDuration, 'alignTrials',p.Results.alignTrials, 'errThr',2);
fprintf('Task struct created: %.1f sec\n', toc)

timing.taskConditionTable = toc


%% 6) Export Task in .mat file (generally for Python)

if p.Results.saveResults
    
    % export Tsync table (for reloading in matlab)
    save(paths.tSyncOutputPath,'Tsync','-v7.3')
    
    % export Task table (for reloading in matlab)
    save(paths.taskTableOutputPathMatlab,'Task','-v7.3')

    % export Task table data (for reading in python)
    exportTaskTable(Task, paths.taskTableOutputPathExport)
    timing.taskDataSaved = toc
end


%% #TODO run functional connectivity analysis 

if p.Results.runConnectivityAnalysis
    warning('runConnectivityAnalysis not implemented')


% % the code below was drawn from run_connected_pairs_analysis_v1 and needs
% adaptation but should be correct 2023-02-13

%     ks.load(loadFeatures=false, loadBatchwise=true);
%     ks.mask_clusters(ks.clusters_good);
%     if numel(ks.clusters_good) < 50
%         warning('only %d good clusters found in dataset: \n%s\nHas manual curation been performed?',numel(ks.clusters_good), ks.pathLeaf)
%         return
%     end
% 
%     % 1) cluster de-duplication preprocessor
%     ntPerBatch = 65536;
%     nBatches = ceil(double(max(ks.spike_times + uint64(100))) / ntPerBatch);
%     ks.batch_sort_order = (1:nBatches)';
%     ks.batch_starts = (uint64(1) : ntPerBatch : (ntPerBatch*uint64(ks.nBatches-1) + uint64(1)))';
% 
%     cda = ClusterDeduplicationAnalysis(ks, ks.clusters_good)
% 
%     delete(gcp('nocreate'))
%     parpool(4); % or some reasonable number of cores, typically I don't use more than 8.
%     cda.detect_duplicates
% 
%     % mask out duplicated clusters
%     keep_clusters = ks.clusters_good(~ismember(ks.clusters_good, cda.remove_duplicate_cluster_ids));
%     ks.mask_clusters(keep_clusters);
% 
%     % 2) run Connected pairs analysis
%     cpa = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', p.Results.jitter_reps);
% 
%     % Compute metrics
%     cpa.computeSmoothedCCGs();
%     cpa.findConnectedPairs;




end


end







