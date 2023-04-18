%% Kinetic/Kinematic Data Class
classdef Kin
    properties (Access = public)
        Fs = 1e3
        data = []
        variableLabels
        alignIndex = 1
        dataXYZ = []
    end
    properties (SetAccess = private)
        nDataPoints
        nVariables
        nTrials
    end
    methods
        % -----------------------------------------------------------------
        % CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = Kin(Fs, data, variableLabels, alignIndex)
            
            if nargin > 0 && isnumeric(Fs) && isscalar(Fs)
                obj.Fs = Fs;
            end
            if nargin > 1 && isnumeric(data) && ndims(data)<=3
                    obj.data = data;
            end
            if nargin > 2 && iscell(variableLabels) && length(variableLabels)==size(data,2)
                obj.variableLabels = variableLabels;
            else
                obj.variableLabels = cellfun(@(x) ['var ' num2str(x)],...
                    num2cell(1:size(obj.data,2)), 'uni', false);
            end
            if nargin > 3 && isnumeric(alignIndex) && isscalar(alignIndex) && alignIndex>0
                obj.alignIndex = alignIndex;
            end
        end
        % -----------------------------------------------------------------
        % SET METHODS
        % -----------------------------------------------------------------
        function obj = set.Fs(obj, Fs)
            if isnumeric(Fs) && isscalar(Fs)
                obj.Fs = Fs;
            else
                error('Invalid sampling frequency assignment')
            end
        end
        function obj = set.data(obj, data)
            if isnumeric(data) && ndims(data)<=3
                obj.data = data;
            else
                error('Invalid data assignment')
            end
        end
        function obj = set.variableLabels(obj, variableLabels)
            if iscell(variableLabels)
                obj.variableLabels = variableLabels;
            else
                error('Invalid variable label assignment')
            end
        end
        function obj = set.alignIndex(obj, alignIndex)
            if isscalar(alignIndex) && alignIndex>0
                obj.alignIndex = alignIndex;
            else
                error('Invalid zero index assignment')
            end
        end
        % -----------------------------------------------------------------
        % GET METHODS
        % -----------------------------------------------------------------
        function varLabel = get.variableLabels(obj)
            nDataVars = size(obj.data,2);
            nLabels = length(obj.variableLabels);
            if nLabels > nDataVars
                varLabel = obj.variableLabels(1:nDataVars);
            elseif nLabels < nDataVars
                varLabel = [obj.variableLabels, cellfun(@(x) ['var ' num2str(x)],...
                    num2cell((1+nLabels):nDataVars), 'uni', false)];
            else
                varLabel = obj.variableLabels;
            end
        end
        function tZerIdx = get.alignIndex(obj)
            tZerIdx = min(obj.alignIndex, size(obj.data,1));
        end
        function nPts = get.nDataPoints(obj)
            nPts = size(obj.data,1);
        end
        function nVar = get.nVariables(obj)
            nVar = size(obj.data,2);
        end
        function nTr = get.nTrials(obj)
            nTr = size(obj.data,3);
        end
        % -----------------------------------------------------------------
        % SELECT DIMENSION
        % -----------------------------------------------------------------
        % time range (dim 1)
        function obj = range(obj,tLim)
            assert(tLim(1)>=0 && tLim(2)<=obj.time(end), 'Time range out of bounds')
            tIdx = (obj.time>=tLim(1) & obj.time<=tLim(2));
            obj.data = obj.data(tIdx,:,:);
            obj.alignIndex = 1-find(tIdx,1);
        end
        % variable (dim 2)
        function obj = var(obj,varIdx)
            assert(all(varIdx<=obj.nVariables),'Variable out of bounds')
            obj.variableLabels = obj.variableLabels(varIdx);
            obj.data = obj.data(:,varIdx,:);
        end
        % trial (dim 3)
        function obj = trial(obj,trNo)
            assert(trNo>=1 && trNo<=obj.nTrials)
            obj.data = obj.data(:,:,trNo);
        end
        % -----------------------------------------------------------------
        % DETECT ONSET
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % PROCESSING
        % -----------------------------------------------------------------
        
        % time vector
        function t = time(obj)
            t = ((1:obj.nDataPoints)-obj.alignIndex)/obj.Fs;
        end
        
        % re-sample
        function obj = resample(obj,newFreq)
            nSteps = size(obj.data,1);
            t = ((1:size(obj.data,1))-obj.alignIndex)'/obj.Fs;
            tNew = ((1:size(obj.data,1))-obj.alignIndex)'/newFreq;
            obj.data = cell2mat(cellfun(@(x) interp1(t,x,tNew),...
                mat2cell(obj.data, nSteps, ones(1,size(obj.data,2)), ones(1,size(obj.data,3))),'uni',false));
            obj.Fs = newFreq;
        end
        
        % trial average
        function [trMean, trSte] = trial_avg(obj)
            trMean = mean(obj.data,3);
            trSte = std(obj.data,[],3)/sqrt(obj.nTrials);
        end
        
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        
        % plot trial
        function plot_trial(obj,varargin)
            P = inputParser;
            addParameter(P,'trial',1:obj.nTrials)
            addParameter(P,'var',1,@isnumeric)
            addParameter(P,'overlay',[],@isnumeric)
            addParameter(P,'clf',false,@islogical)
            addParameter(P,'new',false,@islogical)
            parse(P,varargin{:})
            
            if P.Results.clf
                clf
            end
            if P.Results.new
                figure
            end
            hold on
            
%             shade = linspace(0,0.75,obj.nTrials)'.*ones(obj.nTrials,3);
            
            t = obj.time();
            trNo = P.Results.trial;
            varNo = P.Results.var;
            for ii = 1:length(P.Results.trial)
                plot(t,obj.data(:,varNo,trNo(ii)),'color',0.3*ones(1,3))
            end
            if ~isempty(P.Results.overlay) && length(t)==length(P.Results.overlay)
                plot(t,P.Results.overlay,'r','linewidth',1)
            end
            set(gca,'xlim',t([1 end]))
            xlabel('time (s)')
        end
            
        
        % plot mean
        function plot_mean(obj,varNo,lineSpec)
            if nargin == 2
                lineSpec = repmat({'k-'},1,varNo);
            end
            t = ((1:obj.nDataPoints)-obj.alignIndex)/obj.Fs;
            [trMean,trSte] = trial_avg(obj);
            hold on
            for ii = 1:length(varNo)
                plot(t,trMean(:,varNo(ii)),lineSpec{ii},'linewidth',2)
                patch([t,fliplr(t)], [trMean(:,varNo(ii))'-trSte(:,varNo(ii))',...
                    fliplr(trMean(:,varNo(ii))'+trSte(:,varNo(ii))')], lineSpec{ii}(1),...
                    'EdgeAlpha',0, 'FaceAlpha',0.125)
            end
            set(gca,'xlim',t([1 end]))
            xlabel('time (s)')
            plot(t([1 end]),[0 0],'k--','color',.8*[1 1 1])
%             if length(obj.variableLabels) >= max(varNo)
%                 legend(obj.variableLabels(varNo), 'location', 'best')
%             end
        end
    end
end
        
        