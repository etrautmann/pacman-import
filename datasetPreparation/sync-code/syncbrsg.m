%% SYNCBRSG synchronizes Blackrock with Speedgoat
% This function synchronizes broadband data from Blackrock with Speedgoat,
% trializes the data, and folds it into a table. The broadband data
% vector is only saved in its raw format. The raw data for each trial is
% then filtered, rectified, and downsampled to Speedgoat's time base (1 kHz).
% 
% Synchronization is accomplished by reading the sync signal output from
% Speedgoat. Because this process can incur time reversals/skips in the
% synchronized clock if Blackrock is running fast/slow, the converted (30
% kHz) clock is linearly resampled for each trial. Thus, the time interval
% between sample points will *not* be the same across trials. To conserve
% system memory, only the dt per trial is recorded. To recover the (roughly)
% 30 kHz time vector, simply take t(1):dt:t(end), where t is the time
% vector recorded on Speedgoat at 1 kHz.
%
% This function is computationally intensive. To maximize efficiency,
% variable arguments allow the user to simultaneously synchronize sorted
% spike times without (much) additional overhead. The function also prints
% dynamic updates to allow the user to monitor progess. These updates can
% be suppressed via variable arguments.
%
% SYNTAX
%   newT = syncemg(nsx, T, varargin)
%
% REQUIRED INPUTS
%   nsx (structure) - nsx structure containing data recorded on Blackrock
%   T (table) - "T" table containing continuous acquisition data
%
% OPTIONAL INPUTS: none
%
% VARIABLE INPUTS
%   (...,'fSampSG',<scalar>) - sampling frequency on Speedgoat
%   (...,'emgRaw',<logical>) - if true, folds raw EMG into table
%   (...,'emgFilt',<logical>) - if true, filters raw EMG and saves to table
%   (...,'spikes',<sparse logical>) - array of spike indices of size
%       [timesteps x units]
%   (...,'FiltParams',<struct>) - filter parameters, specified by the
%       following structure fields: order ('ord'), Gaussian SD ('sd'), and 
%       highpass cutoff ('cut')
%   (...,'printStatus',<logical>) - if true, prints function status to
%       command window
%
% OUTPUT
%   newT (table) - new "T" Table with synchronized data
%
%   unitNames (cell array) - names of each sorted unit corresponding to
%       each column in newT. Only valid if syncing spike times.
%
% EXAMPLE:
%   % specify path to xPC data (see LOADBEHAVIOR for additional info)
%   xpcPath = [churchland_data '/speedgoat/2015-09-25/Cousteau_EMG_029'];
%
%   % load continuous acquisition data
%   rawT = loadBehavior(xpcPath);
%
%   % specify path to NS5 file
%   ns5Path = [emgDir '/C_2015_07_30_EMG_BR1_001.ns5'];
%
%   % load NS5
%   ns5 = openNSx(ns5Path);
%
%   % synchronize emg saved on an ns5 file; omit filtered data
%   syncedT = syncemg(ns5.Data, rawT, 'no_filt');
%
% Other m-files required: FILTERGAUSS, STATUSBAR
% Subfunctions: FILTERGAUSS, STATUSBAR
% MAT-files required: none
%
% See also: SYNCSORTED, OPENNSX

% Author: Najja Marshall
% Email: njm2149@cumc.columbia.edu
% Dated: February, 2016

function [newT, Times] = syncbrsg(nsx, T, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'SYNCBRSG';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'nsx', @isstruct)
addRequired(P, 'T', @istable)
addParameter(P, 'dataChan', [], @isnumeric)
addParameter(P, 'syncChan', [], @isnumeric)
addParameter(P, 'saveTag', [], @isscalar)
addParameter(P, 'tol', 1, @isscalar)
addParameter(P, 'data', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'FsSG', 1e3, @isscalar)
addParameter(P, 'emgRaw', false, @islogical)
addParameter(P, 'emgFilt', false, @islogical)
addParameter(P, 'spikes', sparse(0,0), @issparse)
addParameter(P, 'FiltParams', struct('sd',25,'ord',12,'cut',40), @isstruct)
addParameter(P, 'legacy', false, @islogical)
addParameter(P, 'emgPlotParams', struct('plotsOn',false, 'plotDuration', 50,'plotStart',5, 'emgChanLabel',{}), @isstruct) % plots for debuging / inspecting emg

% clear workspace (parser object retains the data while staying small)
parse(P, nsx, T, varargin{:});
clear ans varargin


%% Parse inputs

% Blackrock and Speedgoat sampling frequencies
FsBR = nsx.MetaTags.SamplingFreq;
FsSG = P.Results.FsSG;

% data duration
dataDur = size(nsx.Data,2) / FsBR;

% data and sync channels
syncChan = P.Results.syncChan;
if isempty(syncChan)
    syncChan = nsx.MetaTags.ChannelCount;
end
dataChan = P.Results.dataChan;
if isempty(dataChan)
    dataChan = 1:(syncChan-1);
end

%% Decode sync signal

syncSignal = nsx.Data(syncChan,:) > mean(nsx.Data(syncChan,:));

% number of sample points per ms
nSampPerMs = FsBR/FsSG;

% expected number of samples between edges based on timing code
expectPulseLen = nSampPerMs * [1; 2; 6];
expectPulseLen = [expectPulseLen; nSampPerMs * 106]; % (dropped sync pulses plus inter-pulse gap)

for ii = 1:2
    % find rising and falling edges
    edgeIdx = 1 + find(abs(diff(syncSignal)) > 0.5);
    
    % number of samples between rising and falling edges
    pulseLen = diff(edgeIdx);
    
    if ii == 1
        % remove partial leading or trailing blocks
        syncSignal(1:edgeIdx(find(pulseLen>=expectPulseLen(3),1))) = false;
        syncSignal(edgeIdx(find(pulseLen>=expectPulseLen(3),1,'last')):end) = false;
    end
end

% assign unique pulse lengths to nearest expected value
unqPulseLen = unique(pulseLen');
bins = knnsearch(expectPulseLen, unqPulseLen);

% find the first pulse location in each block
firstPulse = [1, 1+find(ismember(pulseLen, unqPulseLen(bins>2)))];
newBlockIdx = edgeIdx(firstPulse);

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
Block.isCorrupt = Block.isCorrupt | [false;diff(Block.time)>0.2] | Block.time>Block.time(1)+dataDur;

pCorrupt = 100*nnz(Block.isCorrupt)/height(Block);
if pCorrupt > 10
    warning('%.3f%% corrupted sync blocks. Timing estimate may be unreliable',pCorrupt)
end

% clear unnecessary data
% clear syncSignal

%% Convert time base from Blackrock to Speedgoat

Times = struct('blackrock',{nsx.MetaTags.Timestamp(1)+((0:size(nsx.Data,2)-1)/FsBR)'},...
    'speedgoat',{[]});

syncTimes = Times.blackrock(1) + (Block.start(~Block.isCorrupt)/FsBR); % s

histEdges = mean([syncTimes(1:end-1),syncTimes(2:end)],2);
histEdges = [0; histEdges; dataDur]; % s

% assign times to bins
[~, timeBins] = histc(Times.blackrock, histEdges);

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

Times.blackrock(rmvIdx) = [];

% compute the difference between sampled times and the nearest sync signal
nearSyncTimes = syncTimes(timeBins);
dtBlackrock = Times.blackrock - nearSyncTimes;

% infer speedgoat times by adding the differences to speedgoat time stamps
tSG = Block.time(~Block.isCorrupt);
Times.speedgoat = roundn((tSG(timeBins) + dtBlackrock), -6);

% clear unnecessary data
clear nearSyncTimes dtBlackrock timeBins

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
% Blackrock are running simultaneously)
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
    tEdgeIdx = tEdgeIdx - round(0.1 * FsBR);
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
newT.trialTime = cellfun(@(t) t(1):(1/FsBR):t(end), newT.t, 'uni', false);
% newT.trialTime = cellfun(@(tei) Times.blackrock(max(1,tei(1)):max(1,tei(2))), tEdgeIdx, 'uni', false);
% newT.trialTime = cellfun(@(tt) tt-tt(1), newT.trialTime, 'uni', false);

% dt for linearly re-sampled Blackrock times
newT.dt = num2cell(repmat(1/FsBR,height(newT),1));
% newT.dt = cellfun(@(tSg,tIdxBr) diff(tSg)/diff(tIdxBr), tEdges, tEdgeIdx, 'uni', false);


%% Fold EMG into T Table

% remove pre-existing EMG
emgFieldLoc = cellfun(@(s) strncmp(s,'emg',3), fieldnames(newT));
newT(:,emgFieldLoc) = [];

goodIdx = firstGoodTrial:nTrials-1;

if isempty(P.Results.data)
    rawData = nsx.Data(dataChan,:);
    
  
else
    rawData = P.Results.data;
end
clear nsx


% append raw voltage traces (30 kHz, int16)
if P.Results.emgRaw
    
    newT.emgBR = cell(nTrials,1);
    
%     newT.emgBR(goodIdx) = cellfun(@(i) rawData(dataCh, i(1):i(2)), tEdgeIdx(goodIdx), 'uni', false);
    newT.emgBR(goodIdx) = cellfun(@(i,t) rawData(dataChan, i(1)+(0:length(t)-1)), tEdgeIdx(goodIdx), newT.trialTime(goodIdx), 'uni', false);
    
    % vertically align
    newT.emgBR = cellfun(@transpose, newT.emgBR, 'uni', false);
end

% append filtered voltage traces (1 kHz, double)
if P.Results.emgFilt
    
    % extract filter parameters
%     filtOrd = P.Results.FiltParams.ord;
%     cutoff = P.Results.FiltParams.cut;
%     gaussSD = P.Results.FiltParams.sd;
    
    % design highpass filter
%     D = fdesign.highpass('N,F3db', filtOrd, cutoff, FsBR);
%     highPass = design(D, 'butter');
%     
%     newT.emgFilt = cell(nTrials,1);
    
    % evaluate each lead one at a time to reduce memory demands
    for iTr = firstGoodTrial:nTrials
        iTr
        % extract trial indices
        trIdx = tEdgeIdx{iTr}(1):tEdgeIdx{iTr}(2);
        
        % identify appropriate sample points
        sampPts = knnsearch((newT.t{iTr}(1) : newT.dt{iTr} : newT.t{iTr}(end))', makecol(newT.t{iTr}));
        
        % filter EMG 
        lowCut = 40;
        highCut = 500;
        gaussWidthMs = 25;
        [thisTrialEmgFilt, ~] = filterEMG(rawData(:,trIdx), FsBR, 'lowCut', lowCut, 'highCut',highCut, 'gaussWidthMs', gaussWidthMs);
        
        newT.emgFilt{iTr} = thisTrialEmgFilt(:,sampPts);
        
        
%         for jLe = 1:length(dataChan)
%             
%             % high-pass filter
%             emgFilt = filtfilt(highPass.sosMatrix, highPass.ScaleValues, double(rawData(dataChan(jLe), trIdx)));
%             
%             % rectify
%             emgFilt = abs(emgFilt);
%             
%             % low-pass filter with Gaussian
%             emgFilt = filterGaussEmt(emgFilt, 'gaussWidthMs',gaussSD, 'Fs', FsBR);
%             
%             % downsample and append to trial
%             newT.emgFilt{iTr} = [newT.emgFilt{iTr}; emgFilt(sampPts)];            
%         end
        
    end
    
    % vertically align
    newT.emgFilt = cellfun(@transpose, newT.emgFilt, 'uni', false);
end


%% Synchronize sorted spike times

if ~isempty(P.Results.spikes)
    
    % convert spike indices to times
    nUnits = size(P.Results.spikes,2);
    tSpk = cell(1,nUnits);
    for ii = 1:nUnits
        tSpk{ii} = Times.blackrock(P.Results.spikes(:,ii));
    end 
    
    newT.tSpk = cell(height(newT), nUnits);
    
    for iTr = firstGoodTrial:nTrials
        % compute speedgoat to blackrock time gain
        tEdgesBr = Times.blackrock(tEdgeIdx{iTr});
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