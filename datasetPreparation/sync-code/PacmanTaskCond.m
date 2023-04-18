%% Pac-Man Task Condition-Aggregated Data Class
classdef PacmanTaskCond
    properties
        session
        Conditions
        Force
        Emg
        MU
        trialStartIndex
        padDur   % EMT added 2023-04-10
    end
    properties (SetAccess = private)
        unitID
        targetForce
        nConditions
    end
    methods
        
        % -----------------------------------------------------------------
        % CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = PacmanTaskCond(Conditions, padDur)
            
            obj.session = '';
            
            % condition table
            if istable(Conditions)
                obj.Conditions = Conditions;
            else
                error('Input condition table')
            end
            
            % Force data
            obj.Force = repmat(Kin([]), obj.nConditions, 1);
            
            % EMG data
            obj.Emg = repmat(Ephys([]), obj.nConditions, 1);
            
            % Motor unit activity data
            obj.MU = repmat(Neuron([]), obj.nConditions, 1);
            
            % trial start indices
            obj.trialStartIndex = cell(obj.nConditions,1);

            obj.padDur = padDur;
        end
        
        % -----------------------------------------------------------------
        % SET METHODS
        % -----------------------------------------------------------------
        function obj = set.session(obj, sess)
            if ischar(sess) || (iscell(sess) && all(cellfun(@ischar,sess)))
                obj.session = sess;
            else
                error('Invalid session assignment')
            end
        end
        function obj = set.MU(obj, motorUnit)
            if isa(motorUnit,'Neuron')
                obj.MU = motorUnit;
                obj = obj.update_unit_id;
            else
                error('Invalid motor unit assignment')
            end
        end
        % -----------------------------------------------------------------
        % GET METHODS
        % -----------------------------------------------------------------
        function nCond = get.nConditions(obj)
            nCond = size(obj.Conditions,1);
        end
        function targForce = get.targetForce(obj)
            nCond = height(obj.Conditions);
            if nCond > 0
                targForce = cell(nCond,1);
                for ii = 1:nCond
                    dt = 1/obj.Force(ii).Fs;
%                     padDur = (obj.Force(ii).alignIndex-1)/obj.Force(ii).Fs;

                    targForce{ii} = pacmantargfns(obj.Conditions(ii,:),1,obj.padDur,'dt',dt);
                end
            end
            % orient vertically
            targForce = cellfun(@(x) x(:), targForce, 'uni', false);
        end
        % -----------------------------------------------------------------
        % OVERLOAD
        % -----------------------------------------------------------------
        function T = plus(T1,T2)
            
            T = PacmanTaskCond(T1.Conditions);
            
            % THIS ASSUMES EXACTLY MATCHED CONDITIONS!!
            nCond = height(T.Conditions);
            
            % merge force data
            for ii = 1:nCond
                T.Force(ii).Fs = T1.Force(ii).Fs;
                T.Force(ii).data = cat(3,T1.Force(ii).data,T2.Force(ii).data);
                T.Force(ii).variableLabels = T1.Force(ii).variableLabels;
                T.Force(ii).alignIndex = T1.Force(ii).alignIndex;
            end
            
            % merge EMG ??
            
            % merge neural data
            for ii = 1:nCond
                minTrials = min(T1.MU(ii).nTrials,T2.MU(ii).nTrials);
                T.MU(ii).Fs = T1.MU(ii).Fs;
                T.MU(ii).spikes = cat(2,T1.MU(ii).spikes(:,:,1:minTrials),T2.MU(ii).spikes(:,:,1:minTrials));
                T.MU(ii).alignIndex = T1.MU(ii).alignIndex;
%                 T.MU(ii).waveform = cat(3,T1.MU(ii).waveform,T2.MU(ii).waveform); % this would be a good reason to save the normalized EMG
            end
        end
        % -----------------------------------------------------------------
        % KINETICS
        % -----------------------------------------------------------------
        
        % force error
        function err = force_err(obj)
            err = cell(obj.nConditions,1);
            targFn = obj.targetForce;
            for cond = 1:obj.nConditions
                % error function
                targVar = var(targFn{cond});
                if targVar == 0
                    targVar = 1;
                end
                errFcn = @(y) 1 - ((y-targFn{cond})'*(y-targFn{cond}))...
                    /(length(targFn{cond})*targVar);
                % error values
                err{cond} = zeros(obj.Force(cond).nTrials,1);
                for trial = 1:obj.Force(cond).nTrials
                    err{cond}(trial) = errFcn(obj.Force(cond).data(:,1,trial));
                end
            end
        end               
        
        % -----------------------------------------------------------------
        % PLOT TOOLS
        % -----------------------------------------------------------------
        % set figure
        function fh = set_fig(~,code)
            if isempty(code)
                return
            end
            switch code
                case 'clf'
                    clf
                    fh = gcf;
                case 'new'
                    fh = figure;
                otherwise
                    return
            end
        end
        
        % -----------------------------------------------------------------
        % DEFINE RECRUITMENT ORDER
        % -----------------------------------------------------------------
        function [obj, recOrd] = order_mu(obj,rateThresh)
            % slowest increasing ramp condition
            incRampCond = find(cellfun(@(typ,amp) strcmp(typ,'RMP') && amp(1)>0, obj.Conditions.type, obj.Conditions.amplitude));
            [~,maxDur] = max(obj.Conditions.duration(incRampCond));
            slowRamp = incRampCond(maxDur);
            
            % define recruitment order
            refPsth = obj.MU(slowRamp).psth(:,:,1);
            refPsth = mat2cell(refPsth, size(refPsth,1), ones(1,size(refPsth,2)));
            onsetIdx = cellfun(@(x) find(x>rateThresh,1), refPsth, 'uni', false);
            onsetIdx(cellfun(@isempty,onsetIdx)) = {Inf};
            onsetIdx = cell2mat(onsetIdx);
            [~,recOrd] = sort(onsetIdx);
            
            % sort remaining units below threshold by PSTH energy
            psthEner = cellfun(@(x) x'*x, refPsth);
            subThresh = ~isfinite(onsetIdx);
            [~,recOrd2] = sort(psthEner(subThresh),'descend');
            x = recOrd(1+nnz(~subThresh):end);
            recOrd(1+nnz(~subThresh):end) = x(recOrd2);
            
            % re-order spike array and waveform
            for ii = 1:length(obj.MU)
                obj.MU(ii).spikes = obj.MU(ii).spikes(:,recOrd,:);
                if ~isempty(obj.MU(ii).waveform)
                    obj.MU(ii).waveform = obj.MU(ii).waveform(:,:,recOrd);
                end
            end
        end
        % -----------------------------------------------------------------
        % EXPORT TASK DATA
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        function fCond = plot_forces(obj,varargin)
            P = inputParser;
            addParameter(P,'cond',1:obj(1).nConditions,@isnumeric)
            addParameter(P,'fig',[],@(x) ischar(x) && ismember(x,{'clf','new','hold'}))
            addParameter(P,'label',[],@ischar)
            addParameter(P,'subplot',true,@islogical)
            addParameter(P,'trials',false,@islogical)
            
            parse(P,varargin{:})
            
            obj(1).set_fig(P.Results.fig);
            nSess = length(obj);
            nCond = length(P.Results.cond);
            
            if nSess > 1
                cmap = brewermap(nSess,'Set1');
            else
                cmap = [0 0 0];
            end
            
            nCol = ceil(sqrt(nCond));
            nRow = ceil(nCond/nCol);
            
            targFn = obj(1).targetForce;
            fCond = cell(nCond,nSess);
            
            for ii = 1:nCond
                if P.Results.subplot && nCond>1
                    subplot(nRow,nCol,ii)
                else
                    figure
                end
                
                condNo = P.Results.cond(ii);
                
                t = obj(1).Force(condNo).time();
                plot(t,targFn{condNo},'c--','linewidth',2)
                hold on
                plot(t([1 end]),[0 0],'k--','color',.8*[1 1 1]) 
                
                fh = zeros(1,nSess);
                
                for jj = 1:nSess

                    if P.Results.trials
                        frc = squeeze(obj(jj).Force(condNo).data(:,2,:));
                        frc = smooth1D(frc,obj(jj).Force(condNo).Fs,'gau','dim',1);
                        for kk = 1:size(frc,2)
                            fh(jj) = plot(t,frc(:,kk),'k');
                        end
                    else
                        [fMean,fSte] = obj(jj).Force(condNo).trial_avg;
                        fMean = fMean(:,2)';
                        fSte = fSte(:,2)';
                        fMean = smooth1D(fMean,obj(jj).Force(condNo).Fs,'gau','dim',2);
                        fSte = smooth1D(fSte,obj(jj).Force(condNo).Fs,'gau','dim',2);
                        
                        fCond{ii} = fMean;
                        fh(jj) = plot(t,fMean,'color',cmap(jj,:),'linewidth',2);
                        patch([t,fliplr(t)], [fMean-fSte, fliplr(fMean+fSte)], cmap(jj,:), 'EdgeAlpha',0, 'FaceAlpha',0.125)
                    end
                end
                
                if nSess > 1
                    legend(fh,cellfun(@(x) sprintf('session %i',x),num2cell(1:nSess),'uni',false),'location','best')
                end
                
                if ((P.Results.subplot && mod(ii,nCol)==1) || ~P.Results.subplot)
                    ylabel('force (N)')
                end
                xlabel('time (s)')
                title(sprintf('condition %i',condNo))
                
                set(gca,'xlim',t([1 end]))
                set(gca,'ylim',[-2,20])
                pause(0.0001)
                
                box off
            end
            if ~isempty(P.Results.label)
                subplot(nRow,nCol,1+(nRow-1)*nCol)
                th = text(0,0,P.Results.label);
                th.Units = 'normalized';
                th.Position = [-.1 -.2];
                th.FontSize = 12;
            end
        end
        
        % single trial data (Force, EMG, rates)
        function plot_trial(obj,prop,trNo,cond,chan)
            
            if ischar(trNo) && strcmp(trNo,'rand')
                trNo = randsample(1:obj.Force(cond).nTrials,1);
            end

            switch prop
                case 'force'
                    t = ((1:obj.Force(cond).nDataPoints)-obj.Force(cond).tZeroIndex)/obj.Force(cond).Fs;
                    x = obj.Force(cond).data(:,:,trNo);
                    var = 1;
                    labels = [];
                    
                case 'emg'
                    t = ((1:obj.Emg(cond).nDataPoints)-obj.Emg(cond).tZeroIndex)/obj.Emg(cond).Fs;
                    x = obj.Emg(cond).data(:,:,trNo);
                    var = chan;
                    labels = cellfun(@(n) sprintf('channel %i',n), num2cell(var), 'uni', false);
            end
            
            nVar = length(var);
            ax = [];
            for ii = 1:nVar
                if nVar > 1
                    ax(ii) = subplot(nVar,1,ii);
                end
                plot(t,x(:,var(ii)),'k')
                if ii == 1
                    title(sprintf('trial %i',trNo))
                end
                if nVar > 1
                    ylabel(labels{ii})
                end
            end
            if nVar > 1
                linkaxes(ax,'x')
                axis('tight')
            end
        end
        
        % EMG reconstruction
        function plot_emg_recon(obj,varargin)
            P = inputParser;
            addParameter(P,'uoi',[],@(x) isnumeric(x))
            addParameter(P,'trial',[],@(x) isnumeric(x) && isscalar(x))
            addParameter(P,'cond',1,@(x) isnumeric(x) && isscalar(x))
            addParameter(P,'chan',[],@(x) isnumeric(x))
            addParameter(P,'label',[],@ischar)
            addParameter(P, 'simplify', [], @(x) isnumeric(x))
            addParameter(P, 'emph', [], @(x) isnumeric(x))
            parse(P,varargin{:})
            
            % read inputs
            cond = P.Results.cond;
            if isempty(P.Results.chan)
                chan = 1:obj.Emg(cond(1)).nChannels;
            else
                chan = P.Results.chan;
            end
            nChan = length(chan);
            
            % color map
            nUnit = obj.MU(cond).nUnits;
%             cmap = brewermap(nUnit,'Spectral');
            cmap = parula(nUnit);
            cmap = max([0 0 0],cmap-.2);
            
            % unit of interest
            uoi = P.Results.uoi;
            
            if ~isempty(P.Results.emph)
                cmap = [cmap, 0.15*ones(size(cmap,1),1)];
                cmap(P.Results.emph,4) = 1;
                if isempty(uoi)
                    uoi = datasample(P.Results.emph,1);
                end
            end
            
            trNo = P.Results.trial;
            if isempty(trNo)
                if isempty(uoi)
                    trNo = randi(obj.Emg(cond).nTrials);
                else
                    trNo = datasample(find(any(obj.MU(cond).spikes(:,uoi,:),1)), 1);
                end
            end
            
            Fs = obj.Emg(cond).Fs;
            tEMG = obj.Emg(cond).time;

            waveLen = size(obj.MU(cond).waveform,1);
            win = -waveLen/4:waveLen/4-1;
            
            % wave norm
            pkIdx = zeros(nChan,nUnit);
            for ch = 1:nChan
                for un = 1:nUnit
                    wNorm = smooth1D(obj.MU(1).waveform(:,ch,un).^2,obj.MU(1).Fs,'gau','sd',2e-3);
                    [~,pkIdx(ch,un)] = max(wNorm);
                end
            end
            pkIdx(pkIdx < 3/8*waveLen | pkIdx > 5/8*waveLen) = waveLen/2;
%             pkIdx = max(pkIdx,3/8*waveLen);
%             pkIdx = min(pkIdx,5/8*waveLen);
            pkIdx = waveLen/2*ones(nChan,nUnit);
            
            % figure elements
            fh = zeros(1,1+nUnit);
            ax = zeros(1,nChan);
            clf
            
            for ch = 1:nChan
                if nChan > 1
                    ax(ch) = subplot(nChan,1,ch);
                end
                
                X = obj.Emg(cond).data(:,chan(ch),trNo);
                
                if ~isempty(P.Results.simplify)
                    
                    spkLoc = double(findpulses(double(X),Fs,'OutputFormat','logical'));
                    spkZone = smooth1D(spkLoc,Fs,'box','wid',waveLen/(2*Fs)) > 0.5;
                    
                    spkBrd = 1+[find(diff(spkZone(2:end))>0.5 & diff(spkZone(1:end-1))==0),...
                        find(diff(spkZone(2:end))<-0.5 & diff(spkZone(1:end-1))==0)];
                    noiseBrd = [[1;spkBrd(:,2)],[spkBrd(:,1);obj.Emg(cond).nSamples]];
                    
                    spkBrd = mat2cell(spkBrd,ones(1,size(spkBrd,1)),2);
                    noiseBrd = mat2cell(noiseBrd,ones(1,size(noiseBrd,1)),2);
                    
                    sampIdxSpk = cellfun(@(x) x(1):P.Results.simplify(2):(x(2)-1),spkBrd,'uni',false);
                    sampIdxNoise = cellfun(@(x) x(1):P.Results.simplify(1):(x(2)-1),noiseBrd,'uni',false);
                    
                    sampIdx = cell2mat(reshape([sampIdxNoise(1:end-1),sampIdxSpk]',1,length(sampIdxSpk)*2));
                    sampIdx = [sampIdx,sampIdxNoise{end},length(X)];
                else
                    sampIdx = 1:length(X);
                end
                
                fh(1) = plot(tEMG(sampIdx), X(sampIdx),'color',[0 0 0 0.25],'linewidth',2);
                hold on
                
                for un = 1:nUnit
                    wFrm = win + pkIdx(ch,un);
                    sFrm = wFrm - waveLen/2;
                    spkLim = [1-sFrm(1), length(tEMG)-sFrm(end)];
                    
                    spkIdx = find(obj.MU(cond).spikes(:,un,trNo));
                    spkIdx = spkIdx(spkIdx>=spkLim(1) & spkIdx<=spkLim(2));
                    
                    w = obj.MU(cond).waveform(wFrm,chan(ch),un);
                    for kk = 1:length(spkIdx)
                        fh(1+un) = plot(tEMG(sFrm+spkIdx(kk)), w, 'color', cmap(un,:), 'linewidth',2);
                    end
                end
                
                if ch == 1
                    title(sprintf('Condition %i, Trial %i\n\nChannel %i',cond,trNo,chan(ch)))
                else
                    title(sprintf('Channel %i',chan(ch)))
                end
                
                box off
            end
            allLabels = [{'raw'},cellfun(@(n) sprintf('unit %i',n), num2cell(1:nUnit), 'uni', false)];
            legend(fh(fh>0), allLabels(fh>0))
            
            if ch == nChan
                xlabel('time (s)')
                if ~isempty(P.Results.label)
                    th = text(0,0,P.Results.label);
                    th.Units = 'normalized';
                    th.Position = [-.075 -.125];
                    th.FontSize = 12;
                end
            end
            if nChan > 1
                linkaxes(ax,'x')
            end
            set(gca,'xlim',tEMG([1 end]))
        end
        
        % EMG reconstruction
        function plot_emg_resid(obj,varargin)
            P = inputParser;
            addParameter(P,'trial',[],@(x) isnumeric(x) && isscalar(x))
            addParameter(P,'cond',1,@(x) isnumeric(x) && isscalar(x))
            addParameter(P,'chan',1,@(x) isnumeric(x))
            addParameter(P,'label',[],@ischar)
            parse(P,varargin{:})
            
            % read inputs
            cond = P.Results.cond;
            chan = P.Results.chan;
            trNo = P.Results.trial;
            if isempty(trNo)
                trNo = randi(obj.Emg(cond).nTrials);
            end
            
            nUnit = obj.MU(cond).nUnits;
            nChan = length(chan);
            waveLen = size(obj.MU(cond).waveform,1);
            halfWaveLen = waveLen/2;
            win = (0:waveLen-1) - halfWaveLen;
            wEner = squeeze(sum(obj.MU(cond).waveform.^2,1));
            
            tEMG = obj.Emg(cond).time();
            spkBounds = [halfWaveLen, length(tEMG)-halfWaveLen];
            
            % figure elements
            fh = zeros(1,3);
            ax = [];
            clf
            
            for ch = 1:nChan
                if nChan > 1
                    ax(ch) = subplot(nChan,1,ch);
                end
                hold on
                
                y = obj.Emg(cond).data(:,chan(ch),trNo);
                
                fh(1) = plot(tEMG, y, 'color', 0.8*ones(1,3), 'LineWidth', 0.5);
                
                yEst = zeros(size(y));
                
                for un = 1:nUnit
                    if wEner(ch,un) == 0
                        continue
                    end
                    
                    spkIdx = find(obj.MU(cond).spikes(:,un,trNo));
                    spkIdx = spkIdx(spkIdx>spkBounds(1) & spkIdx<spkBounds(2));
                    
                    for kk = 1:length(spkIdx)
                        yEst(win+spkIdx(kk)) = yEst(win+spkIdx(kk)) + obj.MU(cond).waveform(:,ch,un);
                    end
                end
                
                fh(2) = plot(tEMG, yEst, 'b');
                fh(3) = plot(tEMG, double(y)-yEst, 'r');
                
                legend(fh, {'raw','est','resid'})
                
                if ch == 1
                    title(sprintf('Condition %i, Trial %i\n\nChannel %i',cond,trNo,chan(ch)))
                else
                    title(sprintf('Channel %i',chan(ch)))
                end
            end
            
            return
            
            if ch == nChan
                xlabel('time (s)')
                if ~isempty(P.Results.label)
                    th = text(0,0,P.Results.label);
                    th.Units = 'normalized';
                    th.Position = [-.075 -.125];
                    th.FontSize = 12;
                end
            end
            if nChan > 1
                linkaxes(ax,'x')
            end
            set(gca,'xlim',tEMG([1 end]))
        end
    end
    methods (Access = private)
        function obj = update_unit_id(obj)
            obj.unitID = ones(1,size(obj.MU(1).waveform,3));
        end
    end
end