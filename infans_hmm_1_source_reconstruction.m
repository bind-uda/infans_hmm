%% description
% This script reads EEG signals and transfomrs them into 
% the source space by the model provided by Anton Tokariev 
% (https://github.com/babyEEG/AED-exposed-infants).
%
%
% NOTE: Some lines must be changed based on your data path.
% REF : Vidaurre et. al (2018), "Spontaneous cortical activity transiently
%       organises into frequency specific phase-coupling networks".

%% initialization
close all; clear; clc;

subjectList   = setdiff(1:68, [10 19 24:25 45 59:60 66]);
dataPath      = '...';
savePath      = '...';
headModelPath = '...\AED\Head Model\';

% loading source space stuff
load([headModelPath 'Atlas.mat'])
load([headModelPath 'CollapseOperator.mat'])
load([headModelPath 'InverseOperator.mat'])

% leakage correction and dipole resolving initalization
Fs                         = 100;           % sampling rate
[p,q]                      = rat(Fs / 250); % for resampling


%% designing band pass filter (0.4-22 Hz)
lpButter = designfilt('lowpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 22, 'SampleRate', 250, 'DesignMethod', 'butter');     
hpButter = designfilt('highpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 0.4, 'SampleRate', 250, 'DesignMethod', 'butter');

%% transforming the sensor-level signals into the source-level signals
mat_files = cell(length(subjectList),1);
T_all     = cell(length(subjectList),1);

for dataID = 1:length(subjectList)
    tic
    % loads .SET files for each subject  
    EEG_H_AS = pop_loadset([dataPath 'AS_H_' num2str(subjectList(dataID)) '.set']);
    EEG_H_QS = pop_loadset([dataPath 'QS_H_' num2str(subjectList(dataID)) '.set']);
    
    % changes the montage to the common average montage
    % in line with the preprocessing step in
    % (Tokariev et. al, 2019 https://www.nature.com/articles/s41467-019-10467-8)
    EEG_H_AS.data = infans_converts_to_averge_montage(EEG_H_AS.data);
    EEG_H_QS.data = infans_converts_to_averge_montage(EEG_H_QS.data);
    
    % Reorders EEG channels based on the "Inverse Operator" proposed by Anton
    % The order of channels for that model should be as follows:
    % Fp1, Fp2, F7, F3, Fz, F4, F8, T3, C3, Cz, C4, T4, T5, P3, Pz, P4, T6, O1, O2
    eeg_H_AS = infans_reorder_channels(EEG_H_AS.data); % channels X samples
    eeg_H_QS = infans_reorder_channels(EEG_H_QS.data); % channels X samples
    
    clear EEG_H_AS EEG_H_QS

    % removes the DC part from all EEG signals (based on Anton's suggestion)
    eeg_H_AS = eeg_H_AS - mean(eeg_H_AS, 2); % channels X samples
    eeg_H_QS = eeg_H_QS - mean(eeg_H_QS, 2); % channels X samples
    
    % filters, resmaples, and selects 3 minutes
    filtered_H_AS = filtfilt(hpButter,filtfilt(lpButter,eeg_H_AS')); % samples X channels
    filtered_H_QS = filtfilt(hpButter,filtfilt(lpButter,eeg_H_QS')); % samples X channels
    
    filtered_H_AS = resample(filtered_H_AS,p,q); % chnages Fs from 250 to 100 Hz
    filtered_H_QS = resample(filtered_H_QS,p,q); % chnages Fs from 250 to 100 Hz
    
    filtered_H_AS = filtered_H_AS(1: Fs * 180, :); % selects just 3 minutes
    filtered_H_QS = filtered_H_QS(1: Fs * 180, :); % selects just 3 minutes
     
    clear eeg_H_AS eeg_H_QS
    
    % reconstructs source activities
    Source_H_AS = infans_gets_source_signals(filtered_H_AS, InverseOperator, CollapseOperator, Atlas); % samples X parcels
    Source_H_QS = infans_gets_source_signals(filtered_H_QS, InverseOperator, CollapseOperator, Atlas); % samples X parcels
    
    % removes the DC part from all souce activities 
    Source_H_AS = Source_H_AS - mean(Source_H_AS); % samples X parcels
    Source_H_QS = Source_H_QS - mean(Source_H_QS); % samples X parcels

    % concatenates AS and QS for each subject
    Source = vertcat(Source_H_AS, Source_H_QS);

    clear filtered_H_AS filtered_H_QS
     
    mat_files{dataID, 1} = Source;                                    % samples X channels
    T_all{dataID, 1}     = [size(Source_H_AS,1) size(Source_H_QS,1)]; % no. of samples for each epoch
    
    clear Source Source_H_AS Source_H_QS

    elapsed_time = toc;   
    fprintf(['Subject ' num2str(dataID) ' out of ' num2str(length(subjectList)) ', elapsed time: ' num2str(elapsed_time) '\n'])
end
save([savePath 'Data_Source.mat'], 'mat_files', 'T_all');
fprintf('Source reconstruction is done!\n')
fprintf('------------------------------\n')


%% Function Definition
function eegout = infans_converts_to_averge_montage(eegin)
% This function converts the original montage to common average
%
% INPUT: 
%   eegin  : a matrix contains EEG data (channels X samples)
%
% OUTPUT:
%   eegout : the changed version of eegin in term of EEG montage (channels X samples)

    eegout = double(eegin - mean(eegin));
end

function eegout = infans_reorder_channels(eegin)
% This function reorders the channels to make them in line with the
% inverse operator proposed by Anton (see the tutorial file in the following link):
% https://github.com/babyEEG/AED-exposed-infants
% Channels for the inverse operator:
% Fp1, Fp2, F7, F3, Fz, F4, F8, T3, C3, Cz, C4, T4, T5, P3, Pz, P4, T6, O1, O2
%
% The order of channels in the sleep dataset is as follows:
% Fp1, Fp2, F3, F4, C3, C4, P3, P4, O1, O2, F7, F8, T3, T4, T5, T6, Fz, Cz, Pz 
%
%
% INPUT:
%   eegin  : a matrix contains EEG data (channels X samples)
%
% OUTPUT:
%   eegout : the reordered version of eegin (channels X samples)

    inverseOrder = ["Fp1"; "Fp2"; "F7"; "F3"; "Fz"; ...
                    "F4"; "F8"; "T3"; "C3"; "Cz";   ...
                    "C4"; "T4"; "T5"; "P3"; "Pz";   ...
                    "P4"; "T6"; "O1"; "O2"];

    mainOrder    = ["Fp1"; "Fp2"; "F3"; "F4"; "C3"; ...
                    "C4"; "P3"; "P4"; "O1"; "O2"; ...
                    "F7"; "F8"; "T3"; "T4"; "T5"; ...
                    "T6"; "Fz"; "Cz"; "Pz"];

    eegout = zeros(size(eegin));           
    for i = 1:length(inverseOrder)
        index = find(inverseOrder == mainOrder(i));
        eegout(index,:) = eegin(i,:);
    end
    eegout = double(eegout);
end

function Source = infans_gets_source_signals(Sensor, InverseOperator, CollapseOperator, MyAtlas)
% This function reconstructs source activities.
%
% see Tokariev et al., 2019, Cerebral Cortex
%
% INPUTS:
%  Sensor           : filtered EEG data (samples X channels)
%  InverseOperator  : Inverse solution for infant head with N = 19 electrodes
%  CollapseOperator : weights for source signals within 'host' parcels
%  MyAtlas          : assignment of src/verticies to parcels (in MyAtlas.Parcels)
%
% OUTPUTS:
%  Source           : filtered parcel signals (samples X parcels (58))
 
    Np = length(MyAtlas.Parcels);  % number of parcels
    L  = size(Sensor, 1);          % EEG length 
    CollapseOperator = repmat(CollapseOperator, 1, L);

    % source signals  
    src_buf = (InverseOperator * Sensor') .* CollapseOperator;      % (sources X samples)

    % parcel/cortical signals  
    parcel_buf = zeros(Np, L);

    for j = 1:Np
        parcel_buf(j, :) = mean(src_buf(MyAtlas.Parcels{j, 1}, :)); % (parcels X samples) 
    end

    Source = parcel_buf'; % (samples X parcels)   
end