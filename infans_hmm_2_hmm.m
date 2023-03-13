%% description
% This script reads the source activities
% (the output of infans_hmm_2_source.m) and infers the
% HMM-TDE ([1]) models for different numbers of states. The aim of this
% step is to find the optimal number of states [2,3,4]. 
%
% NOTE: Some lines must be changed based on your data path.
% REF : 
%   [1] Vidaurre et. al (2018), "Spontaneous cortical activity transiently
%       organises into frequency specific phase-coupling networks".
%   [2] Baker et. al (2014), "Fast transient networks in spontaneous human brain
%       activity".
%   [3] Stevner et. al (2019), "Discovery of key whole-brain transitions and dynamics during human
%       wakefulness and non-REM sleep".
%   [4] Wang et. al (2020), "Transition and Dynamic Reconfiguration of Whole-Brain Network
%       in Major Depressive Disorder"

%% initialization
% make sure that FIELDTRIP and SMP are NOT in your matlab path
close all; clear; clc;

subjectList   = 1:60;
dataPath      = '...\Source\';
savePath      = '...\Data\Model\';
load([dataPath 'Data_Source.mat'])

% HMM initalization
Fs             = 100; % sampling rate
embeddedlag    = 5;   % embedded lags (in sample) (values = 10[-100ms, 100ms])
PCs            = 0.8; % number of PCs for embedding (to explain about 80% of the variance)
lifeTimeThresh = 3;   % if the length of each state is less than lifeTimeThresh, it will be discarded. (in sample)

% making the options structure for HMM
options                = struct();
options.order          = 0;        % it deactivates MAR model
options.zeromean       = 1;        % the mean of the time series will not be used to drive the states
options.covtype        = 'full';   % to have a full covariance matrix for each state (with off-diagonal elements different from zero)
options.embeddedlags   = -embeddedlag:embeddedlag;
options.pca            = PCs;         
options.Fs             = Fs;
options.onpower        = 0;        % it deactivates estimation of the model based on envelope of signals          
options.standardise    = 1;        % makes the mean of signals zero and the standard deviation one
options.standardise_pc = options.standardise;
options.inittype       = 'HMM-MAR';
options.cyc            = 100;
options.initcyc        = 10;
options.initrep        = 3;
options.verbose        = 1;
options.dropstates     = 0;

% making the options structure for stochastic process
options.BIGNinitbatch      = 0; 
options.BIGNbatch          = 0; 
options.BIGtol             = 1e-7;
options.BIGcyc             = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay           = 1;
options.BIGforgetrate      = 0.7;
options.BIGbase_weights    = 0.9;

% convert type of signals to single
for i = 1:60
    mat_files_single{i,1} = single(mat_files{i,1});
end

%% running the HMM
for K = 3:10
    fprintf(['Working on HMM inference for the state number ' num2str(K) ' ... \n'])
    options.K  = K;
    outputFile = [savePath 'HMM_TDE_K' num2str(K) '_L' num2str(embeddedlag) '_PC' num2str(100 * PCs) '.mat'];
    
    % running HMM-TDE
    tic
    [hmm, Gamma, ~, vpath, ~, ~, ~, ~] = hmmmar(mat_files_single, T_all, options);
    freeEnergy = hmmfe(mat_files_single,T_all,hmm);
    % Gamma = padGamma(Gamma, T_all, options);
    save(outputFile,'hmm','Gamma','vpath', 'freeEnergy', 'options');
    elapsed_time = toc;
    fprintf(['HMM inference is done, elapsed time: ' num2str(elapsed_time) '\n'])

    % computation of temporal features
    % NOTE: modify vpath and Gammma based on different at hand conditions (e.g. AS and QS) and compute each metric for different conditions
    % NOTE: for "LifeTimes" don't change the values related to the length of signals since it is done in the function "getStateLifeTimes"
    fprintf('Working on the calculation of temporal features ... \n')
    metrics = struct();
    
    T = cell2mat(T_all);
    T = T(:);
    T = mat2cell(T, ones(length(T),1));
    Gamma = padGamma(Gamma,T,options);

    % computes temporal metrics for both AS and QS
    metrics.LT = getStateLifeTimes     (vpath, T, options, lifeTimeThresh, 0.66, 0);
    metrics.IT = getStateIntervalTimes (vpath, T, options, lifeTimeThresh, 0.66, 0);
    metrics.FO = getFractionalOccupancy(Gamma, T, options);
    metrics.SR = getSwitchingRate      (Gamma, T, options);
    
    save(outputFile,'metrics','-append')
    
    elapsed_time = toc;
    fprintf(['The calculation of temporal features is done, elapsed time: ' num2str(elapsed_time) '\n'])
    clear hmm Gamma vpath freeEnergy metrics  
end