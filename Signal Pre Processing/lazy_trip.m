% clc;
%clear;
%close all;
% load('task_Day3to5')

restoredefaultpath    % restore default folder for matlab

maindir = pwd;        % keep main path

cd '/Users/niklas/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/FieldTrip' % set up the path of fieldtrip

addpath(pwd)
ft_defaults

cd(maindir)   

% load('task_Day3to5.mat');        % define the data path and its name
dataEEG.label = {'Fp1'; 'F4'; 'F3'; 'FC6'; 'O1'; 'CP6'; 'PO8'; 'PO3'; 'O2'; 'FC1'; 'FC2'; 'CP1'; 'CP2'; 'P3'; 'P4'; 'Fz'; 'Cz'; 'CPz'; 'C3'; 'C4'; 'Fp2'; 'Pz'; 'PO4'};
dataEEG.hdr.label = dataEEG.label;
% Read events
cfg                    = [];
cfg.trialdef.prestim   = 0.1;                   % in seconds
cfg.trialdef.poststim  = 0.2;                   % in seconds
cfg.trialdef.eventtype = 'S12'; % get a list of the available types
% 
% cfg.channel = 'O2';
% epochs = ft_definetrial(cfg);

% cfg.dataset            = data_name;             % set the name of the dataset
% cfg_tr_def             = ft_definetrial(cfg);   % read the list of the specific stimulus

% 
cfg.hpfilter       = 'yes';        % enable high-pass filtering
cfg.lpfilter       = 'yes';        % enable low-pass filtering
cfg.hpfreq         = 1;           % set up the frequency for high-pass filter
cfg.lpfreq         = 40;          % set up the frequency for low-pass filter
cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
cfg.dftfreq        = [50 100 150]; % set up the frequencies for notch filtering
cfg.baselinewindow = [-0.1 -0.02];

% segment data according to the trial definition
    % define the baseline window
data               = ft_preprocessing(cfg, dataEEG);


cfg                    = [];
adc001 = find(strcmp({event.values}, 'S12'));
adc002 = find(strcmp({event.values}, 'S13'));

% get the sample number in the original data
% note that we transpose them to get columns
smp001 = [event(adc001).samples]';
smp002 = [event(adc002).samples]';

pre    =  125;
post   =  250;
offset = -pre; % see ft_definetrial

trl001 = [smp001-pre smp001+post];
trl002 = [smp002-pre smp002+post];

% add the offset
trl001(:,3) = offset;
trl002(:,3) = offset;

trl001(:,4) = 1; % add a column with the condition number
trl002(:,4) = 2; % add a column with the condition number

% concatenate the two conditions and sort them
trl = sortrows([trl001; trl002])


cfg     = [];
cfg.trl = trl;
data_epoch = ft_redefinetrial(cfg,data);
cfg     = [];

ft_rejectvisual(cfg, data_epoch)
%ft_componentanalysis

