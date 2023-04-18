%% PACCOND PacMan Task Conditions for gain switching task
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
%         emt2177@columbia.edu
% Dated:  2021-02-01 

% #TODO: gracefully handle errors when input T table doesn't contain all
% conditions (i.e. in one savetag) EMT 2020-12-16


function [Task, nGlitch, AlignStats, alignIndices] = paccond_gain_switch(T, neuralRecType, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'PACCOND';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'T', @istable);
addRequired(P, 'neuralRecType', @ischar); % specify recording type: {'blackrock','neuropixels','none'}
addParameter(P, 'alignState', 'InTarget', @ischar)
addParameter(P, 'saveTags', [], @isnumeric)
addParameter(P, 'padDur', [], @isnumeric)
addParameter(P, 'filtParams', [], @iscell)
addParameter(P, 'chanNo', [], @isnumeric)
addParameter(P, 'w', [], @isnumeric)
addParameter(P, 'noiseStd',[],@isnumeric)
addParameter(P, 'errThr', 1, @isscalar)
addParameter(P, 'stdThr', 3, @isscalar)
addParameter(P, 'alignTrials', false, @islogical)
addParameter(P, 'FsBr', 30000, @isscalar)   % sampling frequency blackrock
addParameter(P, 'FsNi', 32000, @isscalar)   % sampling frequency National Instruments (I/O for Neuropixels PXIe system) - different from neuropixels headstage
addParameter(P, 'FsNpix',30000, @iscalar)   % sampling freuqnecy for neuropixels headstage
addParameter(P, 'FsSg', 1000, @isscalar)    % sampling frequency speedgoat datalogger


% clear workspace (parser object retains the data while staying small)
parse(P, T, neuralRecType, varargin{:});
clear ans varargin


% Not all of these sampling rates are used on all session types 
FsBr = P.Results.FsBr;  % used for both EMG and some neural recordings
FsNi = P.Results.FsNi;  
FsNpix = P.Results.FsNpix;
FsSg = P.Results.FsSg;
FS_SPIKES_OUT = 1000;


switch neuralRecType
    case 'blackrock'
        FsNeural = FsBr;
    case 'neuropixels'
        FsNeural = FsNpix;
    case 'none'
        FsNeural = [];
    otherwise
        error('neuralRecType must be one of [''blackrock'',''neuropixels'']')
end


%% Constant values

% % condition parameter fields (not used yet)
% COND_PARAMS = {
%     'type'
%     'offset'
%     'amplitude'
%     'duration'
%     'frequency'
%     'power'
%     };

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


% EMT: 2022-07-17 (moving this to loadsession.m)
% %% Add gain from last trial to trialParams
% 
% % gains = zeros(size(T.trialParams));
% T.trialParams{1}.previousGain = nan;
% for ii = 2:height(T) % empy trials cause errors in cellfun, use loop instead
%     try
% %         gains(ii) = T.trialParams{ii}.gain;
%         T.trialParams{ii}.previousGain = T.trialParams{ii-1}.gain;
%     catch 
% %         gains(ii) = nan;
%         T.trialParams{ii}.previousGain = nan;
%     end
% end


%% Initialize task object

% condition fields
trialFields = fieldnames(T.trialParams{1});

% remove extraneous trial fields
trialFields = trialFields(~ismember(trialFields,{'stimDelay'}));

nTrialFields = length(trialFields);
condFields = [{'id'}; trialFields];

% % extract relevant trial parameters
% if ismember('stimElectrode',trialFields)
%     trialParams = cellfun(@(tp) [tp.condNo,tp.stimElectrode,tp.stimCurrent], T.trialParams,'uni',false);
% else
%     trialParams = cellfun(@(tp) tp.condNo, T.trialParams,'uni',false);
% end

% Strip some bad trials
trialMask = ~isnan(T.saveTag);
T = T(trialMask,:);

% extract trial parameters from all trials.
% if perturbations exist in the dataset, add those to the conditionparameter list
if nnz(strcmp(fieldnames(T.trialParams{1}),'pertFlag')) > 0
    trialParams = cellfun(@(tp) [tp.frcPol, tp.type, tp.offset, tp.amplitude, tp.duration, tp.frequency, tp.gain, tp.pertFlag, tp.pertAmp, tp.pertTime], T.trialParams,'uni',false);
else
    trialParams = cellfun(@(tp) [tp.gain, tp.offset, tp.amplitude, tp.duration, tp.frequency], T.trialParams,'uni',false);
end

% trialParams = cellfun(@(tp) [tp.gain], T.trialParams,'uni',false);

% determine number of unique conditions
unqCond = unique(cell2mat(trialParams),'rows');
nCond = size(unqCond,1);
unqCondNo = (1:nCond)';

% map conditions to unique identifier
condID = cellfun(@(tp) unqCondNo(ismember(unqCond,tp,'rows')), trialParams);

if length(unique(condID)) < nCond
    warning('All conditions are not present in the data - nCond: %d, present: %d',nCond, length(unique(condID)))
    nCond = length(unique(condID));
    unqCond = unqCondNo(ismember(unqCondNo,condID));

end

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

% add previous gain into condition table
gain = cellfun(@(tp) tp.gain, T.trialParams);
% previousGain = cellfun(@(tp) tp.previousGain, T.trialParams);

% this was already commented out before removing previousGain fields
% % Condition.previousGain = [];
% % for ii = 1:nCond
% %     condIdx = find(condID == unqCondNo(ii));
% %     Condition.previousGain{ii} =  previousGain(condIdx)
% % end

% Condition.previousGain = [];
% task object
Task = PacmanTaskCond(Condition, P.Results.padDur); %(:,trialFields(ismember(trialFields,COND_PARAMS))));

%% Populate task object

AlignStats = cell(nCond,1);
alignIndices = cell(nCond,1);

maxNRMSE = cell(nCond,1);

pbar = ProgressBar(nCond, 'Processing Condition data', nCond);
for ii = 1:nCond

%     if ii == 20
%         keyboard
%     end

    pbar.update(ii, sprintf('Condition %d of %d',ii,nCond) )
%     pbar.update(ii, 'asdfs' )
    
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
    
    % sample indices within the trial (Speedgoat and Neural)
    trIdxSG = -round(FsSg*padDur(1)) : round(FsSg*(padDur(2)+Condition.duration(ii)));
    trIdxNeural = -round(FsNeural*padDur(1)) : round(FsNeural*(padDur(2)+Condition.duration(ii)));
    
    % check bounds
    % this check makes sure that trials have force data for each of the
    % sample indices in trIdxSG. If too many trials are unbounded, check
    % that padDur(2) isn't too large and beyond the trial boundary 
    isBounded = cellfun(@(ai,fr) diff([ai+trIdxSG(end),length(fr)],[],2)>=0, alignIdx, T.forceRaw(condIdx));
    isBounded = isBounded & cellfun(@(ai,ff) diff([ai+trIdxSG(end),length(ff)],[],2)>=0, alignIdx, T.forceFilt(condIdx));
    if ismember('emgSG',fieldnames(T))
        isBounded = isBounded & cellfun(@(ai,esg) diff([ai+trIdxSG(end),size(esg,2)],[],2)>=0, alignIdx, T.emgSG(condIdx));
    end
    if ismember('emgBR',fieldnames(T))
        isBounded = isBounded & cellfun(@(ai,ebr) diff([ai+trIdxNeural(end),size(ebr,1)],[],2)>=0, alignIdx, T.emgBR(condIdx));
    end
    condIdx = condIdx(isBounded);
    alignIdx = alignIdx(isBounded);
    
    % alignment stats
    AlignStats{ii} = struct('emptyPts',nnz(emptyAlignment), 'unbounded',nnz(~isBounded));
    
    alignIndices{ii} = alignIdx;
    
    colNames = Condition.Properties.VariableNames;
	colInds = ismember(colNames,{'type','gain','blockMembership','offset','amplitude','duration','frequency','decay','power'});
    
    % phase correction for dynamic conditions. 
    % NOTE: this is not appropriate for any analyses requiring detailed
    % timing (i.e. perturbation response, divergence between sets of
    % conditions, etc.)
    if P.Results.alignTrials
        warning('Aligning trials to behavior - impacts timing of neural responses')
        if ~strcmp(Condition.type{ii},'STA')
            MAX_LAG = 0.2;
            padDurTrunc = padDur-2*MAX_LAG;

            targFn = pacmantargfns(Condition(ii,colInds),1,'dt',1/FsSg,'padDur',padDurTrunc);
            tIdx = -round(FsSg*padDurTrunc(1)):round(FsSg*(padDurTrunc(2)+Condition.duration(ii)));
            maxLagSamp = round(MAX_LAG*FsSg);
            lags = -maxLagSamp:maxLagSamp;
            optLag = zeros(length(condIdx),1);
            maxNRMSE{ii} = zeros(length(condIdx),1);
            for trial = 1:length(condIdx)
                normRMSE = zeros(length(lags),1);
                for ll = 1:length(lags)

                    frcAlg = T.forceFilt{condIdx(trial)}(tIdx+alignIdx{trial}+lags(ll))';   
                    if Condition.gain(ii) < 0
                        frcAlg = Condition.frcMax(ii) + Condition.gain(ii)*frcAlg;
                    end
                    normRMSE(ll) = 1 - sqrt(mean((frcAlg-targFn).^2)/var(targFn));
                end
                [~,maxIdx] = max(normRMSE);
                optLag(trial) = lags(maxIdx);
                maxNRMSE{ii}(trial) = normRMSE(maxIdx);
            end

            % shift alignment index
            alignIdx = num2cell(cell2mat(alignIdx)+optLag);
        end
    end
    
    % check improper length trials
    if ismember('emgBR',fieldnames(T))
        badTrials = cellfun(@(x) x*(FsBr/FsSg)+trIdxNeural(end),alignIdx) > cellfun(@(x) size(x,1),T.emgBR(condIdx));
        condIdx(badTrials) = [];
        alignIdx(badTrials) = [];
    end
    
    badTrials = cellfun(@(fr) length(fr),T.forceRaw(condIdx))...
        < cellfun(@(ai) ai+trIdxSG(end),alignIdx);
    condIdx(badTrials) = [];
    alignIdx(badTrials) = [];
    
    % remove inaccurate trials
    targFn = pacmantargfns(Condition(ii,colInds),1,padDur, 'dt',1/FsSg);
    tIdx = -round(FsSg*padDur(1)) : round(FsSg*(padDur(2)+Condition.duration(ii)));
    
    % #TODO: fix this 
%     y = Condition.frcMax(ii) * cell2mat(cellfun(@(ci,ai,lag) T.forceFilt{ci}(tIdx+ai), num2cell(condIdx), alignIdx,'uni',false))';
    
    y = cell2mat(cellfun(@(ci,ai,lag) T.forceFilt{ci}(tIdx+ai), num2cell(condIdx), alignIdx,'uni',false))';
    if Condition.gain(ii) < 0
        y = Condition.frcMax(ii) + Condition.gain(ii)*y;
    end
 
    % 1) filter out trials based on maximum deviation from target
    ii
    absErr = max(abs(y-targFn),[],1)/max(2,max(targFn)-min(targFn));
    badTrials = absErr > P.Results.errThr;
    
    condIdx(badTrials) = [];
    alignIdx(badTrials) = [];
    y(:,badTrials) = [];
    
    % 2) remove imprecise trials
    mu = mean(y,2);
    sd = std(y,[],2);
    badTrials = any(y<(mu-P.Results.stdThr*sd)...
        | y>(mu+P.Results.stdThr*sd),1);
    
    condIdx(badTrials) = [];
    alignIdx(badTrials) = [];


%     Condition.previousGain{ii} =  previousGain(condIdx);
    
    % forces
%     forceRaw = cellfun(@(fr,ai,tp) tp.frcMax*(((MAX_FORCE_POUNDS*NEWTONS_PER_POUND)/tp.frcMax * (fr(ai+trIdxSG)/MAX_FORCE_VOLTS)) - tp.frcOff)',...
%         T.forceRaw(condIdx), alignIdx, T.trialParams(condIdx), 'uni', false);
%     forceFilt = cellfun(@(ff,ai,tp) Condition.frcMax(ii) * ff(ai+trIdxSG)',...
%         T.forceFilt(condIdx), alignIdx, T.trialParams(condIdx), 'uni', false);

    X_FORCE_IND = 1;
    Y_FORCE_IND = 2;
    Z_FORCE_IND = 3;
    
    % leaving these here for backwards compatibility prior to three-axis measurements
    forceRaw = cellfun(@(fr,ai)  fr(Y_FORCE_IND, ai+trIdxSG)',...
        T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
    forceFilt = cellfun(@(ff,ai) ff(Y_FORCE_IND, ai+trIdxSG)',...
        T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);

    ForceRawX = cellfun(@(fr,ai)  fr(X_FORCE_IND, ai+trIdxSG)', T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
    ForceFiltX = cellfun(@(ff,ai) ff(X_FORCE_IND, ai+trIdxSG)', T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);    
    
    ForceRawY = cellfun(@(fr,ai)  fr(Y_FORCE_IND, ai+trIdxSG)', T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
    ForceFiltY = cellfun(@(ff,ai) ff(Y_FORCE_IND, ai+trIdxSG)', T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);    

    ForceRawZ = cellfun(@(fr,ai)  fr(Z_FORCE_IND, ai+trIdxSG)', T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
    ForceFiltZ = cellfun(@(ff,ai) ff(Z_FORCE_IND, ai+trIdxSG)', T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);    

    forces = [{forceFilt}, {forceRaw}];
    forces = cellfun(@(x) permute(cell2mat(x'),[1 3 2]), forces, 'uni', false);
    
    Fxyz = [{ForceFiltX},{ForceFiltY},{ForceFiltZ}];
    Fxyz = cellfun(@(x) permute(cell2mat(x'),[1 3 2]), Fxyz,'uni', false);
    
    Task.Force(ii).data = cell2mat(Fxyz);  % changed from 
    Task.Force(ii).Fs = FsSg;
    Task.Force(ii).alignIndex = 1 + padDur(1)*FsSg;
    Task.Force(ii).variableLabels = {'Fx','Fy','Fz'};
   
    
    % EMG
   
    if ismember('emgFilt',fieldnames(T)) % % already downsampled EMG within syncbrsg
        
        emg = cellfun(@(ebr,ai) ebr(ai+trIdxSG,:), T.emgFilt(condIdx), alignIdx, 'uni', false);
        FsEmg = FsSg;
    
%     if ismember('emgBR',fieldnames(T)) % Blackrock
%         
%         emg = cellfun(@(ebr,ai) ebr(ai*(FsBr/FsSg)+trIdxNeural,:), T.emgBR(condIdx), alignIdx, 'uni', false);
%         fSamp = FsBr;
%         
    elseif ismember('emgSG',fieldnames(T)) % Speedgoat
        
        emg = cellfun(@(esg,ai) esg(:,ai+trIdxSG)', T.emgSG(condIdx), alignIdx, 'uni', false);
        FsEmg = FsSg;   
        
    else
        FsEmg = FsBr;
        emg = [];
    end
    
    if ~isempty(emg)
        Task.Emg(ii).data = cell2mat(permute(emg,[2 3 1]));
        Task.Emg(ii).Fs = FsEmg;
        Task.Emg(ii).alignIndex = 1 + padDur(1)*FsEmg;
    end
    
    

    % spike data
    if ismember('tSpk', fieldnames(T))
        try
            % copy meta data
            Task.MU(ii).Fs = FS_SPIKES_OUT; 
            Task.MU(ii).alignIndex = 1 + padDur(1)*FS_SPIKES_OUT;
            
            % unit count
            nUnits = size(T.tSpk,2);
    
            % align spike times
            % NOTE: FsNI used here because spike times should already be aligned and transformed into NI sample times. 
            % This was the source of the delayed neural activity bug, fixed on 2021-06-06
            tAlign = cellfun(@(tt,ai) tt(ai*(FsNi/FsSg)), T.trialTime(condIdx), alignIdx);  
            spks = cellfun(@(s,tA) s-tA, T.tSpk(condIdx,:), num2cell(repmat(tAlign,1,nUnits)), 'uni', false);
            
            % bound spike times
            tBound = (1+[0 (Task.Force(ii).nDataPoints-1)*(FS_SPIKES_OUT/FsSg)] - Task.MU(ii).alignIndex)/Task.MU(ii).Fs;
            spks = cellfun(@(s) s(s>=tBound(1) & s<=tBound(2)), spks, 'uni', false);
            
            % build sparse spike array
    %         tt = tBound(1):1/Task.MU(ii).Fs:tBound(2);
            tt = tBound(1): 1/1000 : tBound(2); % NOTE: changed this to discritize by 1kHz hardcoded here (instead of FsNeural) for speed
            nTrials = size(spks,1);
            spkTemp = ndSparse.build([1 1 1],false,[length(tt), nUnits, nTrials]);
            
            for un = 1:nUnits
    %             un
                for tr = 1:nTrials
                    spkTemp(:,un,tr) = hist(spks{tr,un},tt);
                end  
            end
            
            Task.MU(ii).spikes = spkTemp;
        catch
            warning("error in processing spikes for condition: %d", ii)
        end
    end
end
pbar.finish('Done')

% swap in the modified version of the conditions table with previousGain
% field added
Task.Conditions = Condition;


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

