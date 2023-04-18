    function [ST] = cousteauNeuropixelDataLUT()
%
%
% EMT 2019-12-10

% Build lookup table of savetags from each day of recording
ST = struct([]);

ST(end+1).date = '2019-11-22';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf 0]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-11-22/raw/pacman-task_c_191122_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-11-25';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-11-25/raw/pacman-task_c_191125_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-11-26';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-11-26/raw/pacman-task_c_191126_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-02';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-02/raw/pacman-task_c_191202_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-03';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-03/raw/pacman-task_c_191203_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-05';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-05/raw/pacman-task_c_191205_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-10';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-10/raw/pacman-task_c_191210_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-11';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-11/raw/pacman-task_c_191211_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-12';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-12/raw/pacman-task_c_191212_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2019-12-13';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-13/raw/pacman-task_c_191213_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

% settling data
ST(end+1).date = '2019-12-13';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-13/raw/pacman-task_c_191213_neu_settling_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';


ST(end+1).date = '2019-12-16';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2019-12-16/raw/pacman-task_c_191216_neu_emg_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';


%%

ST(end+1).date = '2020-01-03';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2020-01-03/raw/pacman-task_c_200103_neu_001.ns6';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';


ST(end+1).date = '2020-01-06';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data2/pacman/Cousteau/2020-01-06/raw/pacman-task_c_200106_neu_emg_001.ns6';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2020-01-08';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data2/pacman/Cousteau/2020-01-08/raw/pacman-task_c_200108_neu_emg_001.ns6';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';



%% diverging conditions task

ST(end+1).date = '2020-01-29';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [50 7100]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2020-01-29/raw/pacman-task_c_200129_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2020-01-31';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2020-01-31/raw/pacman-task_c_200131_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2020-02-04';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 5200]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2020-02-04/raw/pacman-task_c_200204_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2020-03-09';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 5700]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2020-03-09/raw/pacman-task_c_200309_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

ST(end+1).date = '2020-03-11';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'npix128';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '/data1/pacman/Cousteau/2020-03-11/raw/pacman-task_c_200311_neu_001.bin';
ST(end).mapFile = 'neuropixels_primateDemo128_chanMap.mat';

%% Pacman gain switch 

% dates to add:


%% Pacman gain switch perturb


ST(end+1).date = '2021-05-07';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'NHPv1';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '';
ST(end).mapFile = '';
ST(end).saveTags = '0';
ST(end).padDuration = .75; 

ST(end+1).date = '2021-05-10';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'NHPv1';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '';
ST(end).mapFile = '';
ST(end).saveTags = '0';
ST(end).padDuration = .75; 

ST(end+1).date = '2021-05-07';
ST(end).subject = 'Cousteau';
ST(end).probeVersion = 'NHPv1';
ST(end).timeRange = [0 inf]; 
ST(end).rawDataPath = '';
ST(end).mapFile = '';
ST(end).saveTags = '0';
ST(end).padDuration = .75; 

% ST(end+1).date = '';
% ST(end).subject = 'Cousteau';
% ST(end).probeVersion = '';
% ST(end).timeRange = [0 inf]; 
% ST(end).rawDataPath = '';
% ST(end).mapFile = '';
% ST(end).saveTags = '';
% ST(end).padDuration = .75; 




% Date-specific analysis details
ST(end+1).subject = 'cousteau';
ST(end).date = '2021-05-20';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [0];
ST(end).syncChan = 2;
ST(end).syncBit = 8;

ST(end+1).subject = 'cousteau';
ST(end).date = '2021-03-18';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [0];
ST(end).syncChan = 3;
ST(end).syncBit = 8;


ST(end+1).subject = 'cousteau';
ST(end).date = '2021-05-28';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [2];
ST(end).syncChan = 2;
ST(end).syncBit = 8;


ST(end+1).subject = 'cousteau';
ST(end).date = '2020-12-16';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [0];
ST(end).syncChan = 4;
ST(end).syncBit = 2;


ST(end+1).subject = 'cousteau';
ST(end).date = '2021-08-11';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [0:4];
ST(end).syncChan = 3;
ST(end).syncBit = 8;


ST(end+1).subject = 'cousteau';
ST(end).date = '2021-08-14';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [0:4];
ST(end).syncChan = 3;
ST(end).syncBit = 8;



ST(end+1).subject = 'cousteau';
ST(end).date = '2021-08-15';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [0:4];
ST(end).syncChan = 3;
ST(end).syncBit = 8;


ST(end+1).subject = 'cousteau';
ST(end).date = '2022-01-31';
ST(end).gNum = 0;
ST(end).tNum = 0;
ST(end).imecNum = 0;
ST(end).padDuration = .75;
ST(end).SAVE_TAG = [2];
ST(end).syncChan = 3;
ST(end).syncBit = 8;


% Create table
ST = struct2table(ST);

ST.index = [1 : size(ST,1)]';
