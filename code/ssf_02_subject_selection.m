%% Subject and electrode selection for BSEP-rsfMRI analysis
%
% Select stimulation electrodes with available CCEP/BSEP data located in
% Yeo Visual or Control networks and satisfying anatomical criteria.
%
% Author: Maria Guadalupe Yanez Ramos
% Developed with scientific and technical guidance from Dora Hermes
% and the Multimodal Neuroimaging Lab (MNL) team.
% May 2026

clc;
clearvars;

%% Project paths

thisFile = which('ssf_02_subject_selection');

if isempty(thisFile)
    error('Could not locate ssf_02_subject_selection.m on the MATLAB path.');
end

codeDir = fileparts(thisFile);
projectDir = fileparts(codeDir);

localDataPath = fullfile(projectDir, 'data');

addpath(genpath(fullfile(codeDir, 'functions')));
%% Subjects

all_subjects = {'01', '02'};

%% Select subjects and electrodes

vcm_Subjects = ssf_select_subject_electrodes( ...
    localDataPath, ...
    all_subjects);

%% Save selection

outputFile = fullfile( ...
    localDataPath, ...
    'vcm_Subjects.mat');

save(outputFile, 'vcm_Subjects', '-v7.3');

fprintf('\nSaved:\n%s\n', outputFile);