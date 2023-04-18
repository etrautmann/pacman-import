function [msqe] = simulatedneurons(O,C,am,RA,meAn,ntrials,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
clc
P=inputParser;
P.FunctionName = 'simulatedneurons';
addRequired(P,'O',@isstruct);
addRequired(P,'C',@istable);
addRequired(P,'am',@(x) isa(x,'double'))
addRequired(P,'RA',@(x) isa(x,'double'))
addRequired(P,'meAn',@(x) isa(x,'double'))
addRequired(P,'ntrials', @(x) isa(x,'double'))
addOptional(P,'c1',.1,@(x) x>=0)
addParameter(P,'NumberOfLatents',50:50:500,@(x) isa(x,'double'))
addParameter(P,'NumberOfRepetitions',10,@isscalar)
addParameter(P,'SampleSize',0,@isscalar)
addParameter(P,'NumberOfConditions',[3,5,7,9,12],@(x) isa(x,'double'));
addParameter(P,'msqe',cell(1,1,3),@(x) isa(x,'cell'))
addParameter(P,'ProduceGraphs',true,@islogical)
addParameter(P,'CorticalNeurons',true,@islogical)
addParameter(P,'UseSister',false,@islogical)
addParameter(P,'NeuroPixel',false,@islogical)
addParameter(P,'Ordered',false,@islogical)
parse(P,O,C,am,RA,meAn,ntrials,varargin{:})
clear ans varargin
order = [1,2,3,5,7,8,9,4,6,10,11,12];
nreps=P.Results.NumberOfRepetitions;
msqe=P.Results.msqe;
ProduceGraphs= P.Results.ProduceGraphs;
nConditions  = P.Results.NumberOfConditions;
if P.Results.SampleSize == 0
    samplesize = size(meAn,2);
else
    samplesize = P.Results.SampleSize;
end
forCortex = P.Results.CorticalNeurons;
UseSister = P.Results.UseSister;
NPXLOnly = P.Results.NeuroPixel;
ncondix = numel(nConditions);
for ni =1:ncondix
    idrx=O.indices;
    numConditions = nConditions(ni);
    dim=P.Results.NumberOfLatents;
    dim = dim(dim<=samplesize);
%     condindex = find([msqe{:,1,1}] == numConditions,1);
    condindex = find(all(cat(1,msqe{:,1,1}) == repmat([numConditions, Ordered],size(msqe,1),1),2));
    if isempty(condindex)
       oldsize = size(msqe,1);
       condindex = oldsize + 1;
       msqe{condindex,1,1} = [numConditions,P.Results.Ordered];
       msqe{condindex,1,2} = samplesize;
    end
    sampleindex = find([msqe{condindex,:,2}] == samplesize, 1);
    if numConditions ~= 12
        if P.Results.Ordered
            cids = order(1:nConditions(ni));
        else
            cids = randi(12,nConditions(ni));
        end
        gtimes = find(ismember(idrx,cids));
        idrx = idrx(ismember(idrx,cids));
    else
        gtimes = 1:numel(idrx);
    end
    if isempty(sampleindex)
        oldsize = size(msqe,2);
        sampleindex = oldsize + 1;
        msqe{condindex,sampleindex , 1 } = samplesize;
    end
    msor=size(msqe{condindex,sampleindex,3},2); 
    if P.Results.SampleSize ~= 0 || (nConditions(ni) ~= 12)
        if P.Results.SampleSize ~=0
            nmax = size(am,2);
            nperm = randperm(nmax,min(nmax,samplesize));
            Am = am(:,nperm) - meAn;
        else
            nperm = 1:size(am,2);
            Am = am - meAn;
        end
        Am = Am(gtimes,nperm);
        ntrials = ntrials(:,nperm);
        meAn = mean(Am);
        RA = RA(:,nperm);
        if P.Results.CorticalNeurons
            if UseSister
                PsTHE = O.pstheven;
                PsTHO = O.psthodd;
    %             stdpsth = std(Psth);
    %             PsTHE = PsTHE./stdpsth;
    %             PsTHO = PsTHO./stdpsth;
            else
    %             stdpsth = std(PSTH);
                PsTHE=cell2mat(C.even);
    %             PsTHE = PsTHE./stdpsth;
                PsTHO=cell2mat(C.odd);
    %             PsTHO = PsTHO./stdpsth;
            end

        else
    %          stdpsthm = std(PSTHm);
             PsTHE = cell2mat(C_motor.even);
    %          PsTHE = PsTHE./stdpsthm;
             PsTHO=cell2mat(C_motor.odd);
    %          PsTHO = PsTHO./stdpsthm; 
        end
        PsTHE = PsTHE(gtimes,nperm);
        PsTHO = PsTHO(gtimes,nperm);
        pce = pca(PsTHE);
        pco = pca(PsTHO);
        Costco = zeros(samplesize,2);

        for i =1:samplesize
            pje=PsTHE*pce(:,i);
            pjo=PsTHO*pce(:,i);
            pje2=PsTHE*pco(:,i);
            pjo2=PsTHO*pco(:,i);

            Costco(i,1)=corr(pje,pjo);
            Costco(i,2)=corr(pje2,pjo2);
        end
        corac= mean(Costco,2);
    else
        if forCortex
            corac=O.correlation;
        else
            corac = O.correlation_muscle;
        end
    end
    numdim = numel(dim);
    for rep=1:nreps
        dimstorer = zeros(2,numdim);
        for di=1:numdim
            %p=pca(Am);
            clc
            disp(['Trial number: ', num2str(rep), ' for ', num2str(dim(di)),' dimension(s) control using ', num2str(nConditions(ni)), ' condition(s), Ordered: ', num2str(P.Results.Ordered) ])
            p = RandOrthMat(size(Am,2));
            A=(Am)*p(:,1:dim(di))*transpose(p(:,1:dim(di))) + meAn;
            Ahat=zeros(size(A));
            ra=range(A);
            A=A.*(RA./ra);
            A = A - mean(A) + meAn;
            for i=unique(idrx)'
                Ahat(idrx==i,:)=addpoissonnoise(A(idrx==i,:),round(ntrials(i,:)),'shape',1,'sigma',15); %#ok<*PFBNS>
            end
            Ahatnorm=Ahat./(range(Ahat)+2);
            r1=Ahatnorm-mean(Ahatnorm);
%             r1 = r1./std(r1);
            for i=unique(idrx)'
                Ahat(idrx==i,:)=addpoissonnoise(A(idrx==i,:),round(ntrials(i,:)),'shape',1,'sigma',15);
            end
            Ahatnorm=Ahat./(range(Ahat)+2);
            r2=Ahatnorm-mean(Ahatnorm);
%             r2=r2./(std(r2));
            r1 = r1./(max(r1) +5);
            r2 = r2./(max(r2) +5);
            p=pca(r1);
            p2 = pca(r2);
            co=zeros(1,size(r1,2));
            for i=1:size(r1,2)
                try
                   co(i) = 1/2* (corr(r1*p(:,i),r2*p(:,i)) + corr(r1*p2(:,i),r2*p2(:,i)));
                catch
                    disp('caught an error')   
                    disp('resuming operation')
                end
            end
            dimstorer(:,di) = [dim(di),mean((corac-transpose(co)).^2)];
        end
        msqe{condindex,sampleindex,3}(1,msor+1:msor+numdim)=dimstorer(1,:);
        msqe{condindex,sampleindex,3}(2,msor+1:msor+numdim)=dimstorer(2,:);
        clc
        msor = msor + numdim;
        MS=msqe;
        if forCortex
            if UseSister
                save('MSreg','MS')
            elseif NPXLOnly
                save('MSPXL','MS')
            else
                save('MS','MS')
            end
        else
            save('MSM','MS')
        end
    end
end
if ProduceGraphs
    figure('Name','Mean Squared Error Plot')
    err=zeros(1,numel(dim));
    mea=zeros(1,numel(dim));
    Msqe = msqe{condindex,sampleindex,3};
    for i=1:numel(dim)
        err(i)=std(Msqe(2,Msqe(1,:)==dim(i)));
        mea(i)=mean(Msqe(2,Msqe(1,:)==dim(i)));
    end
    errorbar(dim,mea,err,'-o','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
    xlabel('Number of Latent Variables')
    ylabel('Mean Squared Error')
end
if isempty(msqe{1,1,1})
    msqe(1,:,:) = [];
end
MS=msqe;
if forCortex
    if UseSister
        save('MSreg','MS')
    elseif NPXLOnly
        save('MSPXL','MS')
    else
        save('MS','MS')
    end
else
    save('MSM','MS')
end
end
