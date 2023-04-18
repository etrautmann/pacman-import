function [paths] = pacmanPaths(data_root, subject, date, gNum, tNum, imecNum)
% generates paths for neuropixels data files assuming the following folder


% #TODO: update this structure since raw and processed data no longer split
% for local processing EMT 2023-02-14 
%
% structure:
% - data_root
% 	- subject
% 		- raw
% 			- date
% 				- speedgoat
% 				- neuropixels
%                   - recording gNum folder
%       				- probe1 folder
%                       	- file.ap.bin
%                       	- file.lf.bin
%                       	- file.ap.meta
%                       	- file.lf.meta
%                   - nidq.bin
%                  	- nidq.meta
% 				- blackrock
% 		- processed
% 			- date
%               - kilosort-manually-sorted
%                   - recording gNum folder
%                       - probe1 sort folder
%                       - proben sort folder
% 
%
%
% EMT 2021-03-08
% updated 2023-02-14 to include new definitions for
%   paths.tSyncOutputPath, paths.taskTableOutputPathMatlab, paths.taskTableOutputPathExport




% 0) *** prefixes

paths.date = date;
paths.dateShort = date([3,4,6,7,9,10]);

prefix = ['pacman-task_' subject(1) '_' date([3,4,6,7,9,10])];
paths.prefix = prefix;

prefixBehavior = [prefix '_beh'];
paths.prefixBehavior = prefixBehavior;

prefixNpix = ['pacman-task_' subject(1) '_' date([3,4,6,7,9,10]) '_neu'];
paths.prefixNpix = prefixNpix;

% 1) *** behavioral data via speedgoat
sgDataPath = fullfile(data_root, subject, date, 'speedgoat');
warnIfNotExist(sgDataPath)
paths.sgDataPath = sgDataPath;


% 2) *** Neuropixels data files (see folder structure above for documentation of expected paths)]

gNumFolder = [prefixNpix '_g' num2str(gNum)];
probeFolder = [prefixNpix '_g' num2str(gNum) '_imec' num2str(imecNum)];  % #TODO: add support for multiple probes here


% build paths
recordingRoot = fullfile(data_root, subject, date, 'neuropixels', gNumFolder); 
warnIfNotExist(recordingRoot)
paths.recordingRoot = recordingRoot;


% Probe path
nPixProbePath = fullfile(recordingRoot, probeFolder);        % #TODO: add support for multiple probes here
warnIfNotExist(nPixProbePath)

% paths for metrics files
paths.ksMetrics = fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_imec' num2str(imecNum) '_ksMetrics.mat']);
paths.centroidsPath = fullfile(nPixProbePath,[prefixNpix '_g' num2str(gNum) '_imec' num2str(imecNum) '_cluster_centroids.tsv']);

% Raw data paths
npixApPath = fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.ap.bin']);
warnIfNotExist(npixApPath)
paths.npixApPath = npixApPath;

npixApMetaPath = fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.ap.meta']);
warnIfNotExist(npixApMetaPath)
paths.npixApMetaPath = npixApMetaPath;

npixLfpPath  =  fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.lf.bin']);
warnIfNotExist(npixLfpPath)
paths.npixLfpPath = npixLfpPath;

npixLfpMetaPath  =  fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.lf.meta']);
warnIfNotExist(npixLfpMetaPath)
paths.npixLfpMetaPath = npixLfpMetaPath;


% 2.1 NIDAQ I/O
nidaqPath = fullfile(recordingRoot, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.nidq.bin']);
warnIfNotExist(nidaqPath)
paths.nidaqPath = nidaqPath;

nidaqMetaPath = fullfile(recordingRoot, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.nidq.meta']);
warnIfNotExist(nidaqMetaPath)
paths.nidaqMetaPath = nidaqMetaPath;


% Manually-sorted kilosort output 
ksResultsPath = fullfile(data_root, subject, date, 'neuropixels',gNumFolder, probeFolder);
warnIfNotExist(ksResultsPath)
paths.ksResultsPath = ksResultsPath;

% % 3) *** blackrock NSX path
% brDataPath = fullfile(data_root, subject, date, 'blackrock', ['pacman-task_' subject(1) '_' date([3,4,6,7,9,10]) '_emg_001.ns6']);
% warnIfNotExist(brDataPath)
% paths.brDataPath = brDataPath;


% 4) path for intermediate output data
% paths.taskTableOutputPath = fullfile('/Volumes/churchland-locker/eric', subject, 'processed', date, 'mergedTaskData', [prefix '_taskdata.mat']);
paths.tSyncOutputPath = fullfile(recordingRoot, [paths.prefixNpix '_tsync.mat']);
paths.taskTableOutputPathMatlab = fullfile(recordingRoot, [paths.prefixNpix '_tasktable_matlab.mat']);
paths.taskTableOutputPathExport = fullfile(recordingRoot, [paths.prefixNpix '_tasktable.mat']);


% 5) export paths as text file for reading in python code/elsewhere
C = [fieldnames(paths), struct2cell(paths)];
filename = fullfile(paths.recordingRoot,[paths.prefix '_data_paths.xml']);
% writecell(C,filename, 'filetype','text', 'delimiter','\t')
writestruct(paths, filename)

end



function [] = warnIfNotExist(path)
    if exist(path) == 0
        warning([path ' doest not exist'])
    end

end

