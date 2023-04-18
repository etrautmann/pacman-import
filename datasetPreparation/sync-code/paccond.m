%% PACCOND PacMan Task Conditions
% Constructs Pac-Man Task condition object from a trial table.
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
%
% REQUIRED INPUTS
%   reqIn (class) - description
%
% OPTIONAL INPUTS
%   optIn (class) - description
%
% VARIABLE INPUTS
%   (...,'parameterName',value) - description (default: )
%
% OUTPUTS
%   C (table) - condition table
%
% EXAMPLE(S) 
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated:

function [Task, nGlitch, AlignStats, alignIndices] = paccond(T, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'PACCOND';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'T', @istable);
addParameter(P, 'alignState', 'InTarget', @ischar)
addParameter(P, 'saveTags', [], @isnumeric)
addParameter(P, 'padDur', [], @isscalar)
addParameter(P, 'filtParams', [], @iscell)
addParameter(P, 'chanNo', [], @isnumeric)
addParameter(P, 'w', [], @isnumeric)
addParameter(P, 'noiseStd',[],@isnumeric)
addParameter(P, 'errThr', 1, @isscalar)
addParameter(P, 'stdThr', 3, @isscalar)

% clear workspace (parser object retains the data while staying small)
parse(P, T, varargin{:});
clear ans varargin


%% Constant values

% Speedgoat & Blackrock sample rates (Hz)
FS_SG = 1e3;
FS_BR = 3e4;

% condition parameter fields
COND_PARAMS = {
    'type'
    'offset'
    'amplitude'
    'duration'
    'frequency'
    'power'
    };

%% Parse trial table

% write task states to struct
stateNo = T.Properties.UserData.TaskStates(:,1);
stateName = T.Properties.UserData.TaskStates(:,2);
TaskStates = cell2struct(stateNo, stateName);

% verify alignment state exists
assert(ismember(P.Results.alignState, stateName))

% filter by save tags
if ~isempty(P.Results.saveTags)
    T = T(ismember(T.saveTag,P.Results.saveTags), :);
end

% remove glitches
glitched = cellfun(@(ts) ismember(TaskStates.Glitch,ts), T.taskState);
nGlitch = nnz(glitched);
T = T(~glitched,:);


%% Initialize task object

% condition fields
trialFields = fieldnames(T.trialParams{1});

% remove extraneous trial fields
trialFields = trialFields(~ismember(trialFields,{'stimDelay'}));

nTrialFields = length(trialFields);
condFields = [{'id'}; trialFields];

% extract relevant trial parameters
if ismember('stimElectrode',trialFields)
    trialParams = cellfun(@(tp) [tp.condNo,tp.stimElectrode,tp.stimCurrent], T.trialParams,'uni',false);
else
    trialParams = cellfun(@(tp) tp.condNo, T.trialParams,'uni',false);
end

% unique conditions
unqCond = unique(cell2mat(trialParams),'rows');
unqCondNo = (1:length(unqCond))';
nCond = length(unqCondNo);

% map conditions to unique identifier
condID = cellfun(@(tp) unqCondNo(ismember(unqCond,tp,'rows')), trialParams);

% initialize condition table
Condition = cell2table(cell(nCond,1+nTrialFields), 'VariableNames', condFields);

% condition ID
Condition.id = (1:nCond)';

% representative set of trial parameters
Condition(:,trialFields) = cell(nCond,nTrialFields);
for ii = 1:nCond
    trPar = T.trialParams{find(condID == unqCondNo(ii), 1)};
    for jTF = 1:nTrialFields
        Condition.(trialFields{jTF}){ii} = trPar.(trialFields{jTF});
    end
    Condition.type{ii} = char(Condition.type{ii});
    Condition.offset{ii} = Condition.offset{ii} * trPar.frcMax;
    Condition.amplitude{ii} = Condition.amplitude{ii} * trPar.frcMax;
end

% extract scalar entries from cells
for iTF = 1:nTrialFields
    if all(cellfun(@isscalar,Condition.(trialFields{iTF})))
        Condition.(trialFields{iTF}) = cell2mat(Condition.(trialFields{iTF}));
    end
end

% task object
Task = PacmanTaskCond(Condition); %(:,trialFields(ismember(trialFields,COND_PARAMS))));

%% Populate task object

AlignStats = cell(nCond,1);
alignIndices = cell(nCond,1);

% gain settings on FUTEK amplifier
MAX_FORCE_POUNDS = 5;
MAX_FORCE_VOLTS = 5.095;

% unit conversion
NEWTONS_PER_POUND = 4.44822;

maxNRMSE = cell(nCond,1);

for ii = 1:nCond
    ii
    % condition indices
    condIdx = find(condID==unqCondNo(ii));
    
    % alignment point
    if ismember('stimElectrode',trialFields) && (Condition.stim(ii)==1)
        alignIdx = cellfun(@(stm) find(stm==1), T.stim(condIdx),'uni',false);
    else
        alignIdx = cellfun(@(ts) find(ts == TaskStates.(P.Results.alignState),1), T{condIdx,'taskState'}, 'uni', false);
    end
    
    % remove trials without alignment point
    emptyAlignment = cellfun(@isempty,alignIdx);
    condIdx(emptyAlignment) = [];
    alignIdx(emptyAlignment) = [];
    
    % pad duration (can override trial parameter)
    if isempty(P.Results.padDur)
        padDur = Condition.padDur(ii);
    else
        padDur = P.Results.padDur;
    end
    
    % trial indices (Speedgoat and Blackrock)
    trIdxSG = -round(FS_SG*padDur):round(FS_SG*(padDur+Condition.duration(ii)));
    trIdxBR = -round(FS_BR*padDur):round(FS_BR*(padDur+Condition.duration(ii)));
    
    % check bounds
    isBounded = cellfun(@(ai,fr) diff([ai+trIdxSG(end),length(fr)],[],2)>=0, alignIdx, T.forceRaw(condIdx));
    isBounded = isBounded & cellfun(@(ai,ff) diff([ai+trIdxSG(end),length(ff)],[],2)>=0, alignIdx, T.forceFilt(condIdx));
    if ismember('emgSG',fieldnames(T))
        isBounded = isBounded & cellfun(@(ai,esg) diff([ai+trIdxSG(end),size(esg,2)],[],2)>=0, alignIdx, T.emgSG(condIdx));
    end
    if ismember('emgBR',fieldnames(T))
        isBounded = isBounded & cellfun(@(ai,ebr) diff([ai+trIdxBR(end),size(ebr,1)],[],2)>=0, alignIdx, T.emgBR(condIdx));
    end
    condIdx = condIdx(isBounded);
    alignIdx = alignIdx(isBounded);
    
    % alignment stats
    AlignStats{ii} = struct('emptyPts',nnz(emptyAlignment), 'unbounded',nnz(~isBounded));
    
    alignIndices{ii} = alignIdx;
    
    % phase correction for dynamic conditions
    if ~strcmp(Condition.type{ii},'STA')
        MAX_LAG = 0.2;
        padDurTrunc = padDur-2*MAX_LAG;
        targFn = pacmantargfns(Condition(ii,6:11),1,'dt',1/FS_SG,'padDur',padDurTrunc);
        tIdx = -round(FS_SG*padDurTrunc):round(FS_SG*(padDurTrunc+Condition.duration(ii)));
        maxLagSamp = round(MAX_LAG*FS_SG);
        lags = -maxLagSamp:maxLagSamp;
        optLag = zeros(length(condIdx),1);
        maxNRMSE{ii} = zeros(length(condIdx),1);
        for trial = 1:length(condIdx)
            normRMSE = zeros(length(lags),1);
            for ll = 1:length(lags)
                frcAlg = Condition.frcMax(ii) * T.forceFilt{condIdx(trial)}(tIdx+alignIdx{trial}+lags(ll))';
                normRMSE(ll) = 1 - sqrt(mean((frcAlg-targFn).^2)/var(targFn));
            end
            [~,maxIdx] = max(normRMSE);
            optLag(trial) = lags(maxIdx);
            maxNRMSE{ii}(trial) = normRMSE(maxIdx);
        end
        
        % shift alignment index
        alignIdx = num2cell(cell2mat(alignIdx)+optLag);
    end
    
    % check improper length trials
    if ismember('emgBR',fieldnames(T))
        badTrials = cellfun(@(x) x*(FS_BR/FS_SG)+trIdxBR(end),alignIdx) > cellfun(@(x) size(x,1),T.emgBR(condIdx));
        condIdx(badTrials) = [];
        alignIdx(badTrials) = [];
    end
    
    badTrials = cellfun(@(fr) length(fr),T.forceRaw(condIdx))...
        < cellfun(@(ai) ai+trIdxSG(end),alignIdx);
    condIdx(badTrials) = [];
    alignIdx(badTrials) = [];
    
    % remove inaccurate trials
    targFn = pacmantargfns(Condition(ii,6:11),1,'dt',1/FS_SG,'padDur',padDur);
    tIdx = -round(FS_SG*padDur):round(FS_SG*(padDur+Condition.duration(ii)));
    y = Condition.frcMax(ii) * cell2mat(cellfun(@(ci,ai,lag) T.forceFilt{ci}(tIdx+ai), num2cell(condIdx), alignIdx,'uni',false))';
    
    absErr = max(abs(y-targFn),[],1)/max(2,max(targFn)-min(targFn));
    badTrials = absErr > P.Results.errThr;
    
    condIdx(badTrials) = [];
    alignIdx(badTrials) = [];
    y(:,badTrials) = [];
    
    % remove imprecise trials
    mu = mean(y,2);
    sd = std(y,[],2);
    badTrials = any(y<(mu-P.Results.stdThr*sd)...
        | y>(mu+P.Results.stdThr*sd),1);
    
    condIdx(badTrials) = [];
    alignIdx(badTrials) = [];
    
    % forces
    forceRaw = cellfun(@(fr,ai,tp) tp.frcMax*(((MAX_FORCE_POUNDS*NEWTONS_PER_POUND)/tp.frcMax * (fr(ai+trIdxSG)/MAX_FORCE_VOLTS)) - tp.frcOff)',...
        T.forceRaw(condIdx), alignIdx, T.trialParams(condIdx), 'uni', false);
    forceFilt = cellfun(@(ff,ai,tp) Condition.frcMax(ii) * ff(ai+trIdxSG)',...
        T.forceFilt(condIdx), alignIdx, T.trialParams(condIdx), 'uni', false);
    
    forces = [{forceFilt}, {forceRaw}];
    forces = cellfun(@(x) permute(cell2mat(x'),[1 3 2]), forces, 'uni', false);
    
    Task.Force(ii).data = cell2mat(forces);
    Task.Force(ii).Fs = FS_SG;
    Task.Force(ii).alignIndex = 1 + padDur*FS_SG;
    Task.Force(ii).variableLabels = {'filt','raw'};
    
    % EMG
    if ismember('emgBR',fieldnames(T)) % Blackrock
        
        emg = cellfun(@(ebr,ai) ebr(ai*(FS_BR/FS_SG)+trIdxBR,:), T.emgBR(condIdx), alignIdx, 'uni', false);
        fSamp = FS_BR;
        
    elseif ismember('emgSG',fieldnames(T)) % Speedgoat
        
        emg = cellfun(@(esg,ai) esg(:,ai+trIdxSG)', T.emgSG(condIdx), alignIdx, 'uni', false);
        fSamp = FS_SG;   
        
    else
        fSamp = FS_BR;
        emg = [];
    end
    
    if ~isempty(emg)
        Task.Emg(ii).data = cell2mat(permute(emg,[2 3 1]));
        Task.Emg(ii).Fs = fSamp;
        Task.Emg(ii).alignIndex = 1 + padDur*fSamp;
    end
    
    % spike data
    if ismember('tSpk', fieldnames(T))
        
        % copy meta data
        Task.MU(ii).Fs = fSamp; %Task.Emg(ii).Fs;
        Task.MU(ii).alignIndex = 1 + padDur*fSamp;
        
        % unit count
        nUnits = size(T.tSpk,2);

        % align spike times
        tAlign = cellfun(@(tt,ai) tt(ai*(FS_BR/FS_SG)), T.trialTime(condIdx), alignIdx);
        spks = cellfun(@(s,tA) s-tA, T.tSpk(condIdx,:), num2cell(repmat(tAlign,1,nUnits)), 'uni', false);
        
        % bound spike times
        tBound = (1+[0 (Task.Force(ii).nDataPoints-1)*(FS_BR/FS_SG)] - Task.MU(ii).alignIndex)/Task.MU(ii).Fs;
        spks = cellfun(@(s) s(s>=tBound(1) & s<=tBound(2)), spks, 'uni', false);
        
        % build sparse spike array
        tt = tBound(1):1/Task.MU(ii).Fs:tBound(2);
        nTrials = size(spks,1);
        spkTemp = ndSparse.build([1 1 1],false,[length(tt), nUnits, nTrials]);
        
        for un = 1:nUnits
            for tr = 1:nTrials
                spkTemp(:,un,tr) = hist(spks{tr,un},tt);
            end
        end
        
        Task.MU(ii).spikes = spkTemp;
    end
end

% filter
if ~isempty(P.Results.filtParams)
    Task.Emg = Task.Emg.filt(P.Results.filtParams{:});
end

% normalize
if ~isempty(P.Results.noiseStd)
    dataClass = class(Task.Emg(1).data);
    for ii = 1:length(Task.Emg)
        Task.Emg(ii).data = double(Task.Emg(ii).data)./P.Results.noiseStd;
        eval(sprintf('Task.Emg(ii).data = %s(Task.Emg(ii).data);',dataClass))
    end
end

% restrict channels
if ~isempty(P.Results.chanNo)
    for ii = 1:nCond
        Task.Emg(ii) = Task.Emg(ii).chan(P.Results.chanNo);
    end
end

if ~isempty(P.Results.w)
    % add waveforms
    for ii = 1:nCond
        Task.MU(ii).waveform = P.Results.w;
    end
    
    % sort motor units by slow ramp condition
    Task = Task.order_mu(10);
end