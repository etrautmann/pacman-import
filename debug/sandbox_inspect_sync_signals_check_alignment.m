
addpath(genpath('/Users/erictrautmann/Dropbox/shenoy_lab/code/spikes'))

apFile = '/Volumes/emt_ssd_3/data/pacman-task/cousteau/raw/2021-06-03/neuropixels/pacman-task_c_210603_neu_g0/pacman-task_c_210603_neu_g0_imec0/pacman-task_c_210603_neu_g0_t0.imec0.ap.bin'
niFile = '/Volumes/emt_ssd_3/data/pacman-task/cousteau/raw/2021-06-03/neuropixels/pacman-task_c_210603_neu_g0/pacman-task_c_210603_neu_g0_t0.nidq.bin'


[niPath,niFileShort,~] = fileparts(niFile);
metaFile = fullfile(niPath,[niFileShort '.meta']);
niMeta = readSpikeGLXmeta(metaFile)

[apPath,apFileShort,~] = fileparts(apFile);
metaFile = fullfile(apPath,[apFileShort '.meta']);
apMeta = readSpikeGLXmeta(metaFile)


fsap = apMeta.imSampRate
fsni = niMeta.niSampRate

%%


nSamp = 500*30000
fidNi = fopen(niFile)
sni = double(fread(fidNi, [2 nSamp],'*int16'));
sni = sni(2,:);

fidAp = fopen(apFile)
sap = double(fread(fidAp, [385, nSamp],'*int16'));
sap = sap(385,:);
%

% normalize the AP recording of the square wave
sap = sap-min(sap);
sap = sap/max(sap);

% grab only the relevant bit for the NI signal
sni2 = bitget(sni,1,'uint16');

%

figure(1); clf;

plot(sni2,'b')
hold on

plot(sap,'r')


%  Figure out how closely timed the two file start times are


%% find first high going pulse (first known synchronous point

nPulse = numel(find(diff(sni2) > 0));
firstPulseNi = find(diff(sni2) > 0, nPulse)
firstPulseAp = find(diff(sap) > 0, nPulse)



% use simple sample rate conversion backwards to determine file start offsets

expectedFirstPulseAp = firstPulseNi/(fsni/fsap)



figure(2); clf;
plot(expectedFirstPulseAp - firstPulseAp)
hold on
xlabel('seconds into recording')
ylabel('expected-actual # of samples')




%% Check existing file to see if there's drift over the course of the session in single trial spike rates


taskDataFile = '/Users/erictrautmann/data/pacman-task/cousteau/processed/2021-03-18/mergedTaskData/pacman-task_c_210318_taskdata.mat'
% taskDataFile = '/Users/erictrautmann/data/pacman-task/cousteau/processed/2021-05-25/mergedTaskData/pacman-task_c_210525_taskdata.mat'
% taskDataFile = '/Volumes/emt_ssd_3/data/pacman-task/cousteau/processed/2021-05-20/mergedTaskData/pacman-task_c_210520_taskdata.mat'
% taskDataFile = '/Users/erictrautmann/data/pacman-task/cousteau/processed/2020-12-16/mergedTaskData/pacman-task_c_201216_taskdata.mat'

data = load(taskDataFile)





%%


for ii = 3
    
    theseSpikes = full(squeeze(sum(data.spikes{ii},2)));
    theseRates = data.rates{ii};
    theseForces = data.forces{ii};
    targ = data.targetForces{ii};
end



% figure(3); clf;
% plot(theseSpikes(:,1),'g')
% hold on
% plot(theseSpikes(:,floor(size(theseSpikes,2)/2)),'g')
% plot(theseSpikes(:,end),'r')

% figure(4); clf;
% 
% %%
% 
% figure(5); clf;
% subplot(311)
% plot(squeeze(sum(theseRates,2)))
% 
% subplot(312)
% plot(squeeze(theseForces(:,2,:)),'b')
% hold on
% plot(targ,'r','linewidth',4)
% 
% subplot(313)
% imagesc(squeeze(sum(theseRates,2))')
% 
% 
% %%

figure(6); clf;
tiledlayout(3,1)


ax1 = nexttile
plot(squeeze(sum(theseRates,2)),'k')
title('sum of spike rates, all neurons')


ax2 = nexttile
plot(squeeze(theseForces(:,2,:)),'b')
hold on
% plot(targ,'r','linewidth',4)



ax3 = nexttile
imagesc(squeeze(sum(theseRates,2))')
ylabel('trial')
xlabel('time')
title('sum of spike rates, all neurons')



linkaxes([ax1 ax2 ax3],'x')




