%% LOADSESSION
% Loads session data
%
% SYNTAX
%   TT = loadsession(prefix, oldTT, varargin)
%
% REQUIRED INPUTS
%   prefix (string) - data prefix, including file path
%
% OPTIONAL INPUTS
%   oldTT (table) - previously returned trial table. If included, will
%       append only newly added data files
%
% VARIABLE INPUTS
%   (...,'parameterName',value) - description (default: )
%
% OUTPUTS
%   TT (table) - trial table
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

function TT = loadsession(prefix, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'LOADSESSION';

% validation functions
% validFn = @(x) exist('x','var');

% add required, optional, and parameter-value pair arguments
addRequired(P, 'prefix', @ischar)
addOptional(P, 'oldTT', [], @istable)
addParameter(P, 'states', {'TaskState'}, @iscell)

% clear workspace (parser object retains the data while staying small)
parse(P, prefix, varargin{:});
clear ans varargin


%% Load session data

% Load session summary.
fileName = sprintf('%s.summary', prefix);
h = fopen(fileName);
str = fread(h,Inf,'*char')';
fclose(h);
summary = str2struct(str);

states = P.Results.states;
for i = 1:length(states)
    names = fieldnames(summary);
    idx = find(strncmp(states{i},names,length(states{i})));
    fieldName = [states{i},'s'];
    stateNums = cell(length(idx),1);
    stateDefs = cell(length(idx),1);
    for j = 1:length(idx)
        stateNums{j} = str2double(names{idx(j)}(length(states{i})+1:end));
        stateDefs{j} = summary.(names{idx(j)});
    end
    summary = rmfield(summary,names(idx));
    summary.(fieldName) = [stateNums,stateDefs];
end

% Calculate number of files.
fileList = dir([prefix '*']);
lastFile = fileList(end).name;
dotIdx = find(lastFile == '.');
nFiles = str2double(lastFile(dotIdx-4:dotIdx-1));

% % Design velocity/torque filter.
% [filt_b, filt_a] = butter(4, 1/50);

nOldTrials = size(P.Results.oldTT,1);
nTrials = (nFiles-1) - nOldTrials;

if nTrials < 1
    disp('Insufficient data')
    TT = [];
    return
end

% Create array of structs to store data.
R = struct(...
    't', [], ...
    'simTime', [], ...
    'validTrial', [], ...
    'success', [], ...
    'saveTag', [], ...
    'taskState', [], ...
    'forceRaw', [], ...
    'forceFilt', [], ...
    'emgSG', [], ...
    'trialParams', [] ,...
    'reward', [] , ...
    'stim', [] , ...
    'photobox', [] ...
    );
R = repmat(R,1,nTrials); % exclude last file in case it's incomplete

% For each trial...
for trialIdx = 1:nTrials

    % indicator of excluded trials
    excluded = false;

    trialNo = trialIdx + nOldTrials;

    % Load data file.
    fileName = sprintf('%s_%04d.data', prefix, trialNo);
    h = fopen(fileName);
    data = fread(h, Inf, '*uint8');
    fclose(h);
    nTimeBytes = 8;
    nLenBytes  = 2;
    nDataBytes = typecast(uint8(data(nTimeBytes+(1:nLenBytes))),'uint16');
    nBytesPerTrial = nTimeBytes+nLenBytes+nDataBytes;
    data = reshape(data, nBytesPerTrial, []);

    % Extract variables.
    simTime = typecast(reshape(data(1:nTimeBytes,:),1,[]),'double');
    i = nTimeBytes + nLenBytes;

    while i < nBytesPerTrial
        name = [lower(char(data(i+(1:3),1)))','_data'];
        type = char(data(i+4,1));
        len  = double(typecast(uint8(data(i+5:i+6,1)),'uint16'));
        if type == 'D'
            bytesIdx = i+6+(1:len*8);
            eval([name,' = ','typecast(reshape(data(bytesIdx,:),1,[]),''double'');'])
            eval([name,' = reshape(',name,',len,[]);'])
            %             disp(name)
        elseif type == 'U'
            bytesIdx = i+6+(1:len);
            eval([name,' = ','data(bytesIdx,:);'])
        else
            warning('Unrecognized data type.')
        end
        i = bytesIdx(end);
    end

    % Task/motor state.
    R(trialIdx).taskState = tst_data;

    % Check for dropped packets and incomplete trials
    Ts = .001;
    R(trialIdx).validTrial = true;
    %     PreTrial = summary.TaskStates{strcmp(summary.TaskStates(:,2),'PreTrial'),1};
    if ~all(diff(simTime) < 1.5*Ts & diff(simTime) > .5*Ts) % || R(trialIdx).taskState(1) ~= PreTrial
        warning('File %d excluded due to dropped packets.',trialNo)
        excluded = true;
    end
    if ~excluded && R(trialIdx).taskState(end) < 100
        warning('File %d was incomplete and was excluded.',trialNo)
        excluded = true;
    end

    % Load params file.
    if ~excluded
        try
            fileName = sprintf('%s_%04d.params', prefix, trialNo);
            h = fopen(fileName);
            d = fread(h, Inf, '*uint8');
            fclose(h);
            tParamsReceived = typecast(d(1:8),'double');
            if ~ismember(tParamsReceived,simTime)
                warning('File %d was excluded because parameters were logged for the wrong trial.',trialNo)
                excluded = true;
            end
            R(trialIdx).trialParams = str2struct(char(d(9:end))');
        catch
            warning('File %d was excluded due to missing parameters.',trialNo)
            excluded = true;
        end
    end

    % handle excluded trials
    if excluded
        R(trialIdx).validTrial = false;
        R(trialIdx).success = 0;
        R(trialIdx).saveTag = NaN;
        R(trialIdx).trialParams = struct('condNo', 0);
        continue
    end

    % Time.
    R(trialIdx).simTime = simTime;
    R(trialIdx).t = 1000*(R(trialIdx).simTime - R(trialIdx).simTime(1)) + 1; % ms, 1 indexed

    % Raw Force (voltage).
    if exist('for_data','var')
        R(trialIdx).forceRaw = for_data;
    else
        R(trialIdx).forceRaw = fry_data;
    end

    % filter parameters for force filtering
    fc = 25;
    fs = 1000;
    [b,a] = butter(6,fc/(fs/2));

    if exist('frx_data','var')
        R(trialIdx).threeAxisForceRaw = [frx_data; fry_data; frz_data];
        R(trialIdx).threeAxisForceFilt = filtfilt(b,a, (R(trialIdx).threeAxisForceRaw'))';
    else
        R(trialIdx).threeAxisForceRaw = nan(3,size(simTime,2));
    end

    % Filtered Force.
    if exist('fof_data','var')
        R(trialIdx).forceFilt = fof_data;
    else
        R(trialIdx).forceFilt = filtfilt(b,a,fry_data);
    end

    if exist('cur_data','var')
        R(trialIdx).cursorPosition = cur_data;
    else
        R(trialIdx).cursorPosition = nan(size(simTime));
    end


    % EMG. (make cell to ensure all fields have same length)
    d = size(emg_data);
    if nnz(d>1) == 1
        emg_data = {emg_data};
    end
    R(trialIdx).emgSG = emg_data;

    % Reward.
    R(trialIdx).reward = rew_data;

    % Stim pulse.
    if exist('stm_data','var')
        R(trialIdx).stim = stm_data;
    end

    % Photobox.
    R(trialIdx).photobox = frm_data;

    % Result.
    successCode = summary.TaskStates{strcmp(summary.TaskStates(:,2),'Success'),1};
    if R(trialIdx).taskState(end) == successCode
        R(trialIdx).success = 1;
    else
        R(trialIdx).success = 0;
    end

    % Save tag.
    R(trialIdx).saveTag = R(trialIdx).trialParams.saveTag;
    R(trialIdx).trialParams = rmfield(R(trialIdx).trialParams,'saveTag');

    R(trialIdx).trialParams = {R(trialIdx).trialParams};
    
  
    try
        R(trialIdx).previousGain = R(trialIdx-1).trialParams{1}.gain;
    catch 
        R(trialIdx).previousGain = nan;
    end

    try
        R(trialIdx).gain = R(trialIdx).trialParams{1}.gain;
    catch 
        R(trialIdx).gain = nan;
    end



end

% Convert to table and remove invalid trials.
TT = struct2table(R);
TT.Properties.UserData = summary;
fileIndices = 1:nTrials;
if ~isempty(TT)
    validFileNums = nOldTrials + fileIndices(TT.validTrial);
    TT.Properties.UserData.FileNumbers = validFileNums;
    TT = TT(TT.validTrial,:);
end
% TT.validTrial = [];

% Append old table
if ~isempty(P.Results.oldTT)
    %     validFiles = [P.Results.oldTT.Properties.UserData.FileNumbers,...
    %         TT.Properties.UserData.FileNumbers];

    TT = [P.Results.oldTT; TT];
    %     TT.Properties.UserData.FileNumbers = validFiles;
end

end
% ---------------------------------------------------------------------
% str2struct
% ---------------------------------------------------------------------

function s = str2struct(str)

% Mark start and end indices.
defIdx = strfind(str,':=');
semiColIdx = strfind(str,';');
nameStart = [1, semiColIdx(1:end-1)+1];
nameEnd = defIdx-1;
valueStart = defIdx+2;
valueEnd = semiColIdx-1;

% Write name-value pairs to struct.
for i = 1:length(valueStart)
    name = str(nameStart(i):nameEnd(i));
    name(isspace(name)) = [];
    value = str(valueStart(i):valueEnd(i));
    try
        s.(name) = eval(value);
    catch
        warning('There was a problem converting a string to a name-value pair in str2struct. The unknown parameter has been written as a NaN.');
        s.(name) = NaN;
    end
end

end