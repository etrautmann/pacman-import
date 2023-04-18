function [newT, Times] = syncSpeedgoatNidaq(niFile, T, varargin)
% #TODO: documentation below needs updating (after copying from srbrsig.m)
%

%
% SYNTAX
%
% REQUIRED INPUTS
%
% OPTIONAL INPUTS: none
%
% VARIABLE INPUTS

%
% OUTPUT
%
% EXAMPLE:

% Other m-files required: FILTERGAUSS, STATUSBAR
% Subfunctions: FILTERGAUSS, STATUSBAR
% MAT-files required: none
%
% See also: SYNCSORTED, OPENNSX

% Author: Eric Trautmann, modified from srbrsig by Najja Marshall
% the only major changes from Najja's code are to enable the function to
% import from NIDAQ files for use with Neuropixels recordsings
% (instead of NSx code for Blackrock recordings)


%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'SYNCSGNIDAQ';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'niFile', @ischar)
addRequired(P, 'T', @istable)
addParameter(P, 'dataChan', [], @isnumeric)
addParameter(P, 'SGsyncChan', 3, @isnumeric)  % channel number in .nidq.bin file used for the speedgoat synchronization signal
addParameter(P, 'SGsyncBit',8,@isnumeric)
addParameter(P, 'saveTag', [], @isscalar)
addParameter(P, 'tol', 2, @isscalar)
addParameter(P, 'data', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'FsSG', 1000, @isscalar)
addParameter(P, 'FsImec', 30000, @isscalar)
addParameter(P, 'emgRaw', false, @islogical)
addParameter(P, 'emgFilt', false, @islogical)
addParameter(P, 'spikes', sparse(0,0), @issparse)
addParameter(P, 'FiltParams', struct('sd',25,'ord',12,'cut',40), @isstruct)
addParameter(P, 'legacy', false, @islogical)

% clear workspace (parser object retains the data while staying small)
parse(P, niFile, T, varargin{:});
clear ans varargin


%% Parse inputs


% import meta file for NIDAQ
[niPath,niFileShort,~] = fileparts(niFile);
metaFile = fullfile(niPath,[niFileShort '.meta']);
if exist(metaFile,'file') == 2
    niMeta = readSpikeGLXmeta(metaFile);
else
    error('No meta file for NIDAQ file: %s', niFileShort)
end


FsNI = niMeta.sRateHz;
FsSG = P.Results.FsSG;

% data duration
dataDur = niMeta.fileTimeSecs; %size(nsx.Data,2) / FsNI;

% data and sync channels
% syncChan = P.Results.syncChan;
% if isempty(syncChan)
%     syncChan = nsx.MetaTags.ChannelCount;
% end
% dataChan = P.Results.dataChan;
% if isempty(dataChan)
%     dataChan = 1:(syncChan-1);
% end



% extract sync signal 
sprintf('Extracting sync signal from NI file')
tic
digitalInputWord = extractSyncChannelFromFile(niFile, niMeta.nChans, P.Results.SGsyncChan);
syncSignal = bitget(digitalInputWord, P.Results.SGsyncBit);

clear digitalInputWord
sprintf('Done extracting sync signal from NI file')
toc

%% Deal with edge cases

% The first two days of recording used a ridiculously high sampling rate, can downsample this
if FsNI == 1000000
    DOWNSAMPLE_FACTOR = 100;
%     syncDat = downsample(syncDat,DOWNSAMPLE_FACTOR);
    FsNI = FsNI/DOWNSAMPLE_FACTOR;
end

% Decode sync signal
% syncThresh = mean([max(syncDat(1:100000)), min(syncDat(1:100000))]);
% syncSignal = syncDat > syncThresh; 

%%

% number of sample points per ms
nSampPerMsNI = FsNI/FsSG;

% expected number of samples between edges based on timing code
expectPulseLen = nSampPerMsNI * [1; 2; 6];
expectPulseLen = [expectPulseLen; nSampPerMsNI * 106]; % (dropped sync pulses plus inter-pulse gap)

for ii = 1:2
    % find rising and falling edges
    edgeIdx = 1 + find(abs(diff(syncSignal)) > 0.5);
    
    % number of samples between rising and falling edges
    pulseLen = diff(edgeIdx);
    
%     if ii == 1
%         % remove partial leading or trailing blocks
%         syncSignal(1:edgeIdx(find(pulseLen>=expectPulseLen(3)+TOLERANCE,1))) = false;
%         syncSignal(edgeIdx(find(pulseLen>=expectPulseLen(3)+TOLERANCE,1,'last')):end) = false;
%     end
end

% assign unique pulse lengths to nearest expected value
unqPulseLen = unique(pulseLen');
bins = knnsearch(expectPulseLen, unqPulseLen);

% find the first pulse location in each block
firstPulse = [1, 1+find(ismember(pulseLen, unqPulseLen(bins>2)))];
newBlockIdx = edgeIdx(firstPulse);

% eliminate last block, which most likely is incomplete
newBlockIdx(end) = [];

% aggregate code into timing blocks
nBlocks = length(newBlockIdx);
Block = repmat(struct('start',[],'code',[],'time',NaN,'isCorrupt',false),nBlocks,1);
powers = 2.^(0:31);
for ii = 1:nBlocks
    Block(ii).start = newBlockIdx(ii);
    Block(ii).code = pulseLen((firstPulse(ii)-1)+(1:2:64));
    
    % infer time stamp from uncorrupted blocks
    if any(min(abs(expectPulseLen(1:2)-Block(ii).code)) > P.Results.tol)
        Block(ii).isCorrupt = true;
    else
        binCode = round((Block(ii).code-expectPulseLen(1))/expectPulseLen(1));
        Block(ii).time = (binCode*powers')/10;
    end
end
Block = struct2table(Block);

% mark as corrupted any blocks with unusually large jumps in time
Block.isCorrupt = Block.isCorrupt | [false; diff(Block.time)>0.2 ] | [diff(Block.time)<0; false] | Block.time>Block.time(1)+dataDur;

pCorrupt = 100*nnz(Block.isCorrupt)/height(Block);
if pCorrupt > 10
    warning('%.3f%% corrupted sync blocks. Timing estimate may be unreliable',pCorrupt)
end

% clear unnecessary data
% clear syncSignal

%% Convert time base from NIDAQ to Speedgoat

% Times = struct('blackrock',{nsx.MetaTags.Timestamp(1)+((0:size(nsx.Data,2)-1)/FsNI)'},...
%     'speedgoat',{[]});

nSamp = niMeta.fileSizeBytes/2/niMeta.nSavedChans;

% handle edge case for high sample rate day early in testing
if niMeta.sRateHz == 1000000
    Times = struct('nidaq',{((0:DOWNSAMPLE_FACTOR:nSamp)/niMeta.sRateHz)'}, 'speedgoat', {[]});
else
    Times = struct('nidaq',{((0:nSamp)/FsNI)'}, 'speedgoat', {[]});
end

syncTimes = Times.nidaq(1) + (Block.start(~Block.isCorrupt)/FsNI); 

histEdges = mean([syncTimes(1:end-1),syncTimes(2:end)],2);
histEdges = [0; histEdges; dataDur]; % s

% assign times to bins
[~, timeBins] = histc(Times.nidaq, histEdges);

% remove data outside the bounds of the sync signal using a forwards and
% backwards sweep
rmvIdx = [];
ii = 1;
while timeBins(ii) < 1
    rmvIdx = [rmvIdx, ii];
    ii = ii + 1;
end
ii = length(timeBins);
lenSyncTimes = length(syncTimes);
while timeBins(ii) > lenSyncTimes
    rmvIdx = [rmvIdx, ii];
    ii = ii - 1;
end

timeBins(rmvIdx) = [];

Times.nidaq(rmvIdx) = [];

% compute the difference between sampled times and the nearest sync signal
nearSyncTimes = syncTimes(timeBins);
dtNidaq = Times.nidaq - nearSyncTimes;

% infer speedgoat times by adding the differences to speedgoat time stamps
tSG = Block.time(~Block.isCorrupt);
Times.speedgoat = roundn((tSG(timeBins) + dtNidaq), -6);

% clear unnecessary data
clear nearSyncTimes dtNidaq timeBins

%% Process behavioral data table

% trial count
nTrials = height(T);

% filter by successful and valid trials and save tag
goodTrials = T.validTrial;
if ismember('success',fieldnames(T))
    goodTrials = goodTrials & T.success;
    T.success = [];
end
if ~isempty(P.Results.saveTag)
    goodTrials = goodTrials & T.saveTag==P.Results.saveTag;
end
T = T(goodTrials,:);

% clear unnecessary data
clear goodTrials

%% Identify trial start and end indices

% update trial count
nTrials = height(T);

% trial start and end times
if ismember('simTime', fieldnames(T))
    tEdges = [cellfun(@(t) t(1), T.simTime), cellfun(@(t) t(end), T.simTime)];
else
    tEdges = [cellfun(@(t) t(1), T.t), cellfun(@(t) t(end), T.t)];
end
tEdges = roundn(tEdges, -3);

% find the first "good" trial (i.e. the first trial when Speedgoat and
% nidaq are running simultaneously)
ii_tEd = 1;
jj_tEd = 1;
while tEdges(ii_tEd, jj_tEd) < Times.speedgoat(1)
    jj_tEd = mod(jj_tEd + 1, 2);
    
    if jj_tEd == 0
        jj_tEd = 2;
    else
        ii_tEd = ii_tEd + 1;
    end
end

firstGoodTrial = ii_tEd + 1;

% match trial start and end indices ("Ed" = "edge"; "Sg" = Speedgoat)
ii_tSg = 1;
ii_tEd = firstGoodTrial;
jj_tEd = 1;

nTimeSteps = length(Times.speedgoat);

tEdgeIdx = zeros(size(tEdges));

while ii_tSg <= nTimeSteps && ii_tEd <= nTrials
    
    if Times.speedgoat(ii_tSg) >= tEdges(ii_tEd,jj_tEd)
        tEdgeIdx(ii_tEd,jj_tEd) = ii_tSg;
        
        jj_tEd = mod(jj_tEd+1, 2);
        if jj_tEd == 0
            jj_tEd = 2;
        else
            ii_tEd = ii_tEd + 1;
        end
    end
    ii_tSg = ii_tSg + 1;
end

if ~P.Results.legacy
    tEdgeIdx = tEdgeIdx - round(0.1 * FsNI);
end

% convert data to cell
tEdgeIdx = mat2cell(tEdgeIdx, ones(1,nTrials), 2);

% build new table
newT = T;
clear T

% convert speedgoat times to seconds and initialize to 0
% newT.t = cellfun(@(t) round((t-t(1))/fSampSG, ceil(log10(fSampSG))), newT.t, 'uni', false);
newT.t = cellfun(@(t) round(t/FsSG, ceil(log10(FsSG))), newT.t, 'uni', false);

% trial times
newT.trialTime = cellfun(@(t) t(1):(1/FsNI):t(end), newT.t, 'uni', false);
% newT.trialTime = cellfun(@(tei) Times.nidaq(max(1,tei(1)):max(1,tei(2))), tEdgeIdx, 'uni', false);
% newT.trialTime = cellfun(@(tt) tt-tt(1), newT.trialTime, 'uni', false);

% dt for linearly re-sampled Blackrock times
newT.dt = num2cell(repmat(1/FsNI,height(newT),1));
% newT.dt = cellfun(@(tSg,tIdxBr) diff(tSg)/diff(tIdxBr), tEdges, tEdgeIdx, 'uni', false);


%% Fold EMG into T Table

if ~isempty(P.Results.data)
	error('LOADING EMG from NIDAQ NOT IMPLEMENTED YET, need to load .nidaq.bin file and parse EMG if desired')
end

% EMG import and filtering code is not yet implemented as of 2020-12-14 (EMT)

% % remove pre-existing EMG
% emgFieldLoc = cellfun(@(s) strncmp(s,'emg',3), fieldnames(newT));
% newT(:,emgFieldLoc) = [];
% 
% goodIdx = firstGoodTrial:nTrials-1;
% 
% if isempty(P.Results.data)
% %     rawData = nsx.Data;
% else
%     rawData = P.Results.data;
% end
% clear nsx
% 
% % append raw voltage traces (30 kHz, int16)
% if P.Results.emgRaw
%     
%     newT.emgBR = cell(nTrials,1);
%     
% %     newT.emgBR(goodIdx) = cellfun(@(i) rawData(dataCh, i(1):i(2)), tEdgeIdx(goodIdx), 'uni', false);
%     newT.emgBR(goodIdx) = cellfun(@(i,t) rawData(dataChan, i(1)+(0:length(t)-1)), tEdgeIdx(goodIdx), newT.trialTime(goodIdx), 'uni', false);
%     
%     % vertically align
%     newT.emgBR = cellfun(@transpose, newT.emgBR, 'uni', false);
% end
% 
% % append filtered voltage traces (1 kHz, double)
% if P.Results.emgFilt
%     
%     % extract filter parameters
%     filtOrd = P.Results.FiltParams.ord;
%     cutoff = P.Results.FiltParams.cut;
%     gaussSD = P.Results.FiltParams.sd;
%     
%     % design highpass filter
%     D = fdesign.highpass('N,F3db', filtOrd, cutoff, FsNI);
%     highPass = design(D, 'butter');
%     
%     newT.emgFilt = cell(nTrials,1);
%     
%     % evaluate each lead one at a time to reduce memory demands
%     for iTr = firstGoodTrial:nTrials
%         
%         % extract trial indices
%         trIdx = tEdgeIdx{iTr}(1):tEdgeIdx{iTr}(2);
%         
%         % identify appropriate sample points
%         sampPts = knnsearch((newT.t{iTr}(1):newT.dt{iTr}:newT.t{iTr}(end))', newT.t{iTr});
%         
%         for jLe = 1:length(dataChan)
%             
%             % high-pass filter
%             emgFilt = filtfilt(highPass.sosMatrix, highPass.ScaleValues, double(rawData(dataChan(jLe), trIdx)));
%             
%             % rectify
%             emgFilt = abs(emgFilt);
%             
%             % low-pass filter with Gaussian
%             emgFilt = filterGauss(emgFilt, gaussSD * FsNI/FsSG);
%             
%             % downsample and append to trial
%             newT.emgFilt{iTr} = [newT.emgFilt{iTr}; emgFilt(sampPts)];            
%         end
%     end
%     
%     % vertically align
%     newT.emgFilt = cellfun(@transpose, newT.emgFilt, 'uni', false);
% end


%% Synchronize sorted spike times

if ~isempty(P.Results.spikes)
    
    % convert spike indices to times
    nUnits = size(P.Results.spikes,2);
    tSpk = cell(1,nUnits);
    for ii = 1:nUnits
        lastInd = min(size(Times.nidaq,1), size(P.Results.spikes,1));
        tSpk{ii} = Times.nidaq(P.Results.spikes(1:lastInd,ii));  %index into times of each sample, using sample numbers of each spike
    end 
    
    newT.tSpk = cell(height(newT), nUnits);
    
    for iTr = firstGoodTrial:nTrials
        % compute speedgoat to blackrock time gain
        tEdgesBr = Times.nidaq(tEdgeIdx{iTr});
        gain = 1; %diff(tEdges{iTr})/(diff(tEdgesBr));
        
        for jUn = 1:nUnits
            % find spikes that fall within trial bounds
            trialSpikes = tSpk{jUn}(tSpk{jUn} >= tEdgesBr(1) & tSpk{jUn} <= tEdgesBr(2));
            
            % record spike times in T table
            newT.tSpk{iTr,jUn} = gain * (trialSpikes - tEdgesBr(1));
        end
    end
end


%% Post-processing

% remove trials preceding the start of the sync signal
if firstGoodTrial > 1
    newT(1:firstGoodTrial-1,:) = [];
end