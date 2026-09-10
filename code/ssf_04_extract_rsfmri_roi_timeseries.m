%% ssf_04_extract_rsfmri_roi_timeseries
%
% Create gray-matter-constrained electrode ROIs in resting-state fMRI
% space and extract rs-fMRI time series for individual electrode contacts
% and bipolar stimulation pairs.
%
% INPUTS
%   data/vcm_Subjects.mat
%   derivatives/electrodes2rsfMRI
%   derivatives/tedana
%   derivatives/fmriprep
%
% REQUIRED PREPROCESSING
%   The fMRIPrep ribbon mask must already be transformed to tedana space
%   and binarized:
%
%   *_desc-ribbon_mask_tedanaSpace_bin.nii.gz
%
% OUTPUTS
%   derivatives/maskROIS/sub-<label>
%   derivatives/timeSeries/sub-<label>/timeSeries.mat
%
% The saved MAT file contains:
%   timeSeries     - voxelwise rs-fMRI signals for individual contacts
%   timeSeriesStim - voxelwise rs-fMRI signals for stimulation-pair ROIs
%
% Author: Maria Guadalupe Yanez Ramos
% Developed with scientific and technical guidance from Dora Hermes
% and the Multimodal Neuroimaging Lab (MNL) team.
% May 2026


clearvars;
clc;
tic;


%% ------------------------------------------------------------------------
% Paths and setup
% -------------------------------------------------------------------------

scriptPath = which('ssf_04_extract_rsfmri_roi_timeseries');

if isempty(scriptPath)
    error( ...
        ['Could not determine the location of ' ...
         'ssf_04_extract_rsfmri_roi_timeseries.m.']);
end

codeDir = fileparts(scriptPath);
projectDir = fileparts(codeDir);

localDataPath = fullfile(projectDir, 'data');
derivativesPath = fullfile(localDataPath, 'derivatives');

addpath(genpath(fullfile(codeDir, 'functions')));


%% ------------------------------------------------------------------------
% Load stimulation-pair selection from Stage 02
% -------------------------------------------------------------------------

subjectsFile = fullfile( ...
    localDataPath, ...
    'vcm_Subjects.mat');

if ~isfile(subjectsFile)
    error( ...
        ['vcm_Subjects.mat was not found. ' ...
         'Run ssf_02_subject_selection first.']);
end

tmp = load(subjectsFile, 'vcm_Subjects');
vcm_Subjects = tmp.vcm_Subjects;

all_subjects = fieldnames(vcm_Subjects)';


%% ------------------------------------------------------------------------
% Analysis settings
% -------------------------------------------------------------------------

ses_label = 'compact3T01';
ieeg_ses_label = 'ieeg01';

rest_run = '01';

circleRadius = 4;


%% ------------------------------------------------------------------------
% Process subjects
% -------------------------------------------------------------------------

for iSubj = 1:numel(all_subjects)

    subj = all_subjects{iSubj};

    fprintf( ...
        '\nProcessing sub-%s (%d of %d)\n', ...
        subj, ...
        iSubj, ...
        numel(all_subjects));


    %% --------------------------------------------------------------------
    % Input paths
    % ---------------------------------------------------------------------

    sesDir = ['ses-' ses_label];

    %% Tedana resting-state fMRI

    tedanaDir = fullfile( ...
        derivativesPath, ...
        'tedana', ...
        ['sub-' subj]);

    fmri_name = fullfile( ...
        tedanaDir, ...
        ['sub-' subj ...
         '_' sesDir ...
         '_task-rest_run-' rest_run ...
         '_desc-optcomAccepted_bold.nii.gz']);

    if ~isfile(fmri_name)
        error( ...
            'Tedana BOLD file not found for sub-%s:\n%s', ...
            subj, ...
            fmri_name);
    end

    fmri = niftiRead(fmri_name);

    %% Electrode coordinates in rs-fMRI T1 space

    elecName = fullfile( ...
        derivativesPath, ...
        'electrodes2rsfMRI', ...
        ['sub-' subj], ...
        ['sub-' subj ...
         '_ses-' ieeg_ses_label ...
         '_space-rsT1_electrodes.tsv']);

    if ~isfile(elecName)
        error( ...
            'Electrode file not found for sub-%s:\n%s', ...
            subj, ...
            elecName);
    end

    elec_table_rs = readtable( ...
        elecName, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'TreatAsEmpty', {'N/A', 'n/a'});


    %% Gray-matter ribbon mask in tedana space

    ribbonMask = fullfile( ...
        derivativesPath, ...
        'fmriprep', ...
        ['sub-' subj], ...
        sesDir, ...
        'anat', ...
        ['sub-' subj '_' sesDir ...
         '_desc-ribbon_mask_tedanaSpace_bin.nii.gz']);

    if ~isfile(ribbonMask)
        error( ...
            ['Gray-matter ribbon mask in tedana space was not found ' ...
             'for sub-%s:\n%s'], ...
            subj, ...
            ribbonMask);
    end

    fmriRibbon = niftiRead(ribbonMask);

    nGrayMatterVoxels = ...
        nnz(fmriRibbon.data == 1);

    fprintf( ...
        'Gray-matter voxels: %d\n', ...
        nGrayMatterVoxels);


    %% --------------------------------------------------------------------
    % Check fMRI and ribbon dimensions
    % ---------------------------------------------------------------------

    fmriSpatialSize = ...
        size(fmri.data);

    ribbonSize = ...
        size(fmriRibbon.data);

    if ~isequal( ...
            fmriSpatialSize(1:3), ...
            ribbonSize(1:3))

        error( ...
            ['Tedana BOLD and ribbon mask dimensions do not match ' ...
             'for sub-%s.'], ...
            subj);
    end


    %% --------------------------------------------------------------------
    % Create gray-matter-constrained ROIs around individual contacts
    % ---------------------------------------------------------------------

    roiDir = fullfile( ...
        derivativesPath, ...
        'maskROIS', ...
        ['sub-' subj]);

    if ~exist(roiDir, 'dir')
        mkdir(roiDir);
    end

    el2analyze = cell(0, 1);

    nElectrodes = height(elec_table_rs);

    for i = 1:nElectrodes

        fprintf( ...
            'Electrode ROI %d of %d\n', ...
            i, ...
            nElectrodes);

        thisChName = ...
            char(string(elec_table_rs.name(i)));

        electrodePosition = [ ...
            elec_table_rs.x(i), ...
            elec_table_rs.y(i), ...
            elec_table_rs.z(i)];

        [niROI, ~, ~] = ...
        ssf_position2resliced_rs_image( ...
            electrodePosition, ...
            fmri_name, ...
            circleRadius);

        %% Binarize spherical ROI

        niROI.data(niROI.data >= 0.5) = 1;
        niROI.data(niROI.data < 0.5) = 0;


        %% Restrict ROI to gray matter

        niROI.data = ...
            niROI.data .* ...
            double(fmriRibbon.data == 1);


        %% Save only ROIs containing gray-matter voxels

        saveName = fullfile( ...
            roiDir, ...
            ['ROI_' thisChName '.nii']);

        niftiWrite( ...
            niROI, ...
            saveName);

        if nnz(niROI.data) > 0

            el2analyze{end + 1, 1} = ...
                thisChName; %#ok<SAGROW>

        else

            fprintf( ...
                '  No gray-matter voxels for %s.\n', ...
                thisChName);

        end

        clear niROI

    end

    fprintf( ...
        'Electrode ROIs with gray-matter voxels: %d of %d\n', ...
        numel(el2analyze), ...
        nElectrodes);


    %% --------------------------------------------------------------------
    % Create ROI for each bipolar stimulation pair
    % ---------------------------------------------------------------------

    subjectInfo = ...
        vcm_Subjects.(subj);

    nStimPairs = ...
        height(subjectInfo);

    for i = 1:nStimPairs

        el1 = ...
            char(string(subjectInfo.stimEl1(i)));

        el2 = ...
            char(string(subjectInfo.stimEl2(i)));

        stimPair = ...
            [el1 '-' el2];

        roiEl1File = fullfile( ...
            roiDir, ...
            ['ROI_' el1 '.nii']);

        roiEl2File = fullfile( ...
            roiDir, ...
            ['ROI_' el2 '.nii']);


        % Preserve the historical requirement that both stimulation
        % contacts have gray-matter-constrained ROIs.
        if ~isfile(roiEl1File) || ~isfile(roiEl2File)

            error( ...
                ['Individual electrode ROI missing for stimulation ' ...
                 'pair %s in sub-%s.'], ...
                stimPair, ...
                subj);

        end

        ROIel1 = ...
            niftiRead(roiEl1File);

        ROIel2 = ...
            niftiRead(roiEl2File);

        stimPairROI = ...
            ROIel1;

        stimPairROI.data = ...
            ROIel1.data + ROIel2.data;

        stimPairROI.data( ...
            stimPairROI.data > 1) = 1;


        % Require at least one gray-matter voxel in the stimulation-pair ROI

        if nnz(stimPairROI.data) == 0

            error( ...
                ['No gray-matter voxels found for stimulation pair %s ' ...
                 'in sub-%s.'], ...
                stimPair, ...
                subj);

        end


        saveName = fullfile( ...
            roiDir, ...
            ['ROI_' stimPair '.nii']);

        niftiWrite( ...
            stimPairROI, ...
            saveName);

        clear ROIel1 ROIel2 stimPairROI
    end

    %% --------------------------------------------------------------------
    % Reshape rs-fMRI data
    % ---------------------------------------------------------------------

    fmriSize = ...
        size(fmri.data);

    if numel(fmriSize) ~= 4
        error( ...
            'Expected 4-D rs-fMRI data for sub-%s.', ...
            subj);
    end

    nVoxels = ...
        prod(fmriSize(1:3));

    nVolumes = ...
        fmriSize(4);

    fmri_data = reshape( ...
        fmri.data, ...
        nVoxels, ...
        nVolumes);


    %% --------------------------------------------------------------------
    % Extract time series from individual-electrode ROIs
    % ---------------------------------------------------------------------

    timeSeries = struct( ...
        'name', {}, ...
        'ts', {});

    for i = 1:numel(el2analyze)

        elName = ...
            el2analyze{i};

        openName = fullfile( ...
            roiDir, ...
            ['ROI_' elName '.nii']);

        ni = ...
            niftiRead(openName);

        roi_data = reshape( ...
            ni.data, ...
            nVoxels, ...
            1);

        timeSeries(i).name = ...
            elName;

        % Rows = ROI voxels
        % Columns = fMRI volumes
        timeSeries(i).ts = ...
            fmri_data(roi_data > 0, :);

        clear ni roi_data

    end


    %% --------------------------------------------------------------------
    % Extract time series from stimulation-pair ROIs
    % ---------------------------------------------------------------------

    timeSeriesStim = struct( ...
        'name', {}, ...
        'ts', {});

    for i = 1:nStimPairs

        el1 = ...
            char(string(subjectInfo.stimEl1(i)));

        el2 = ...
            char(string(subjectInfo.stimEl2(i)));

        stimPair = ...
            [el1 '-' el2];

        openName = fullfile( ...
            roiDir, ...
            ['ROI_' stimPair '.nii']);

        ni = ...
            niftiRead(openName);

        roi_data = reshape( ...
            ni.data, ...
            nVoxels, ...
            1);

        timeSeriesStim(i).name = ...
            string(stimPair);

        % Rows = ROI voxels
        % Columns = fMRI volumes
        timeSeriesStim(i).ts = ...
            fmri_data(roi_data > 0, :);

        clear ni roi_data

    end


    %% --------------------------------------------------------------------
    % Save time series
    % ---------------------------------------------------------------------

    timeSeriesDir = fullfile( ...
        derivativesPath, ...
        'timeSeries', ...
        ['sub-' subj]);

    if ~exist(timeSeriesDir, 'dir')
        mkdir(timeSeriesDir);
    end

    saveName = fullfile( ...
        timeSeriesDir, ...
        'timeSeries.mat');

    save( ...
        saveName, ...
        'timeSeries', ...
        'timeSeriesStim', ...
        '-v7.3');

    fprintf( ...
        'Saved rs-fMRI ROI time series for sub-%s.\n', ...
        subj);

end


fprintf('\n');
toc;