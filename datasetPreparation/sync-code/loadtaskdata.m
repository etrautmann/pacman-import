%% Constants
DATE = '2019-12-16';
DATA_PATH = '/Volumes/Churchland-labshare/Jumanji/pacman-task/cousteau/raw/';
NSX_FILE_NAME = 'pacman-task_c_191216_neu_emg_001.ns6';
SAVE_TAG = 1;

%% load behavior data

filePrefix = ['pacman-task_c_' DATE([3,4,6,7,9,10]) '_'];
T = loadsession([DATA_PATH, DATE '/speedgoat/' filePrefix 'beh']);

%% load nsx data

nsx = openNSX([DATA_PATH, DATE '/blackrock/' NSX_FILE_NAME]);

%% sync kilosort results

tSpk = readNPY([dataPath 'spike_times.npy']);
clus = readNPY([dataPath 'spike_clusters.npy']);
label = tdfread([dataPath 'cluster_groups.csv'],'\t');

label.group = mat2cell(label.group, ones(size(label.group,1),1), size(label.group,2));
label.group = cellfun(@(x) strtrim(x), label.group, 'uni',false);

goodUnit = cellfun(@(x) strcmp(x,'good'), label.group);
multiUnit = cellfun(@(x) strcmp(x,'mua'), label.group);

spkIdx = cellfun(@(x) tSpk(clus==x),...
    num2cell(label.cluster_id(goodUnit | multiUnit)), 'uni',false);
spkIdx = cellfun(@unique, spkIdx, 'uni', false);

s = cell2mat(cellfun(@(idx) sparse(double(idx),1,true,size(nsx.Data,2),1), spkIdx', 'uni',false));

Tsync = syncbrsg(nsx, T, 'spikes',s);

%% load task data

Task = paccond(Tsync,'saveTags',SAVE_TAG,'padDur',1);