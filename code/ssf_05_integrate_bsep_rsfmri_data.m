%% ssf_05_integrate_bsep_rsfmri_data
%
% Integrate BSEP/CRP results with electrode localization information and
% resting-state fMRI ROI time series.
%
% INPUTS
%   data/vcm_Subjects.mat
%   derivatives/whimstim/sub-<label>/*_preproc_CRP.mat
%   derivatives/loc_info/sub-<label>/loc_info.mat
%   derivatives/timeSeries/sub-<label>/timeSeries.mat
%
% OUTPUTS
%   derivatives/vcm_data/sub-<label>/<subject>_<stim-pair>.mat
%   updated data/vcm_Subjects.mat
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

scriptPath = which('ssf_05_integrate_bsep_rsfmri_data');

if isempty(scriptPath)
    error( ...
        ['Could not determine the location of ' ...
         'ssf_05_integrate_bsep_rsfmri_data.m.']);
end

codeDir = fileparts(scriptPath);
projectDir = fileparts(codeDir);

localDataPath = fullfile(projectDir, 'data');
derivativesPath = fullfile(localDataPath, 'derivatives');

addpath(genpath(fullfile(codeDir, 'functions')));


%% ------------------------------------------------------------------------
% Load subject/stimulation-pair information
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

z_thresh = 5;

% Historical preprocessing rule for the stimulation-pair rs-fMRI
% time series.
nDiscardVolumes = 10;


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
    % Load rs-fMRI ROI time series
    % ---------------------------------------------------------------------

    timeSeriesFile = fullfile( ...
        derivativesPath, ...
        'timeSeries', ...
        ['sub-' subj], ...
        'timeSeries.mat');

    if ~isfile(timeSeriesFile)
        error( ...
            'timeSeries.mat not found for sub-%s:\n%s', ...
            subj, ...
            timeSeriesFile);
    end

    tmpTs = load( ...
        timeSeriesFile, ...
        'timeSeries', ...
        'timeSeriesStim');

    timeSeries = tmpTs.timeSeries;
    timeSeriesStim = tmpTs.timeSeriesStim;


    %% --------------------------------------------------------------------
    % Load loc_info
    % ---------------------------------------------------------------------

    locInfoFile = fullfile( ...
        derivativesPath, ...
        'loc_info', ...
        ['sub-' subj], ...
        'loc_info.mat');

    if ~isfile(locInfoFile)
        error( ...
            'loc_info.mat not found for sub-%s:\n%s', ...
            subj, ...
            locInfoFile);
    end

    tmpLoc = load(locInfoFile, 'loc_info');
    loc_info = tmpLoc.loc_info;


    %% --------------------------------------------------------------------
    % Verify stimulation-pair correspondence
    % ---------------------------------------------------------------------

    subjectInfo = vcm_Subjects.(subj);

    nStimPairs = height(subjectInfo);

    expectedPairNames = ...
        string(subjectInfo.stimEl1) + "-" + ...
        string(subjectInfo.stimEl2);

    actualPairNames = ...
        string({timeSeriesStim.name})';

    if numel(actualPairNames) ~= nStimPairs || ...
            ~isequal(expectedPairNames(:), actualPairNames(:))

        error( ...
            ['Stimulation-pair names in timeSeriesStim do not match ' ...
             'vcm_Subjects for sub-%s.'], ...
            subj);
    end


    %% --------------------------------------------------------------------
    % Initialize stimulation-pair mean rs-fMRI time series
    % ---------------------------------------------------------------------

    vcm_Subjects.(subj).meanTimeSeries = ...
        cell(nStimPairs, 1);


    %% --------------------------------------------------------------------
    % Process stimulation pairs
    % ---------------------------------------------------------------------

    for j = 1:nStimPairs

        el1 = ...
            char(string(subjectInfo.stimEl1(j)));

        el2 = ...
            char(string(subjectInfo.stimEl2(j)));

        stim_pair = ...
            [el1 '-' el2];

        fprintf( ...
            'Stim pair %d of %d: %s\n', ...
            j, ...
            nStimPairs, ...
            stim_pair);


        %% ----------------------------------------------------------------
        % Load CRP results from Stage 03
        % -----------------------------------------------------------------

        crpFile = fullfile( ...
            derivativesPath, ...
            'whimstim', ...
            ['sub-' subj], ...
            [subj '_' stim_pair '_preproc_CRP.mat']);

        if ~isfile(crpFile)
            error( ...
                'CRP file not found for sub-%s, pair %s:\n%s', ...
                subj, ...
                stim_pair, ...
                crpFile);
        end

        tmpCrp = load(crpFile, 'crp');
        vcm_data = tmpCrp.crp;


        %% ----------------------------------------------------------------
        % Mean rs-fMRI time series for stimulation contacts
        % -----------------------------------------------------------------

        timeSeriesNames = string({timeSeries.name});

        idxEl1 = find( ...
            timeSeriesNames == string(el1), ...
            1);

        idxEl2 = find( ...
            timeSeriesNames == string(el2), ...
            1);


        % Preserve the historical calculation: concatenate voxelwise
        % time series from both stimulation-contact ROIs before averaging.
        %
        % If only one contact contains gray-matter voxels, use the available
        % contact. This extends the historical calculation to stimulation pairs
        % retained because at least one contact satisfies the selection criteria.

        if ~isempty(idxEl1) && ~isempty(idxEl2)

            tsStim = [ ...
                timeSeries(idxEl1).ts; ...
                timeSeries(idxEl2).ts];

        elseif ~isempty(idxEl1)

            tsStim = ...
                timeSeries(idxEl1).ts;

            fprintf( ...
                '  Only %s contributes gray-matter voxels to %s.\n', ...
                el1, ...
                stim_pair);

        elseif ~isempty(idxEl2)

            tsStim = ...
                timeSeries(idxEl2).ts;

            fprintf( ...
                '  Only %s contributes gray-matter voxels to %s.\n', ...
                el2, ...
                stim_pair);

        else

            error( ...
                ['Neither stimulation contact has gray-matter rs-fMRI ' ...
                 'voxels for pair %s in sub-%s.'], ...
                stim_pair, ...
                subj);

        end


        tsStimmean = ...
            mean(tsStim, 1);

        if size(tsStimmean, 2) <= nDiscardVolumes
            error( ...
                ['Not enough rs-fMRI volumes after discarding the ' ...
                 'first %d volumes for sub-%s, pair %s.'], ...
                nDiscardVolumes, ...
                subj, ...
                stim_pair);
        end

        vcm_Subjects.(subj).meanTimeSeries{j} = ...
            tsStimmean(nDiscardVolumes + 1:end);

        %% ----------------------------------------------------------------
        % Add localization, CRP summary, and rs-fMRI information
        % -----------------------------------------------------------------

        for i = 1:numel(vcm_data)

            recCh = ...
                char(string(vcm_data(i).use_channels_name));


            %% loc_info data

            recCh_li = find( ...
                ismember(string(loc_info.name), string(recCh)), ...
                1);

            if ~isempty(recCh_li)

                vcm_data(i).x = ...
                    loc_info.x(recCh_li);

                vcm_data(i).y = ...
                    loc_info.y(recCh_li);

                vcm_data(i).z = ...
                    loc_info.z(recCh_li);

                vcm_data(i).hemisphere = ...
                    loc_info.hemisphere(recCh_li);

                vcm_data(i).Destrieux_label = ...
                    loc_info.Destrieux_label(recCh_li);

                vcm_data(i).Destrieux_label_text = ...
                    loc_info.Destrieux_label_text(recCh_li);

                vcm_data(i).yeo7VisualControl = ...
                    loc_info.yeo7VisualControl(recCh_li);

                vcm_data(i).gm_wm_relativeDistance = ...
                    loc_info.gm_wm_relativeDistance(recCh_li);

            end


            %% CRP-derived response metrics

            cod = ...
                vcm_data(i).cod;

            med_cod = ...
                median(cod, 'omitnan');

            mad_cod = ...
                mad(cod, 1);

            z_robust = ...
                abs(cod - med_cod) ./ ...
                (1.4826 * mad_cod);

            cod_clean = ...
                cod(z_robust <= z_thresh);

            vcm_data(i).median_cod = ...
                median(cod, 'omitnan');

            vcm_data(i).mean_cod_noOutliers = ...
                mean(cod_clean, 'omitnan');

            vcm_data(i).numtrials = ...
                numel(cod_clean);

            vcm_data(i).Vsnr = ...
                mean(vcm_data(i).crp_parms.Vsnr);

            vcm_data(i).tR = ...
                vcm_data(i).crp_parms.tR;

            vcm_data(i).p_val = ...
                vcm_data(i).crp_projs.p_value_tR;


            %% Significant BSEP response criterion

            if vcm_data(i).p_val <= 0.05 && ...
                    vcm_data(i).tR > 0.020 && ...
                    vcm_data(i).Vsnr >= 1 && ...
                    vcm_data(i).mean_cod_noOutliers >= 0.1

                vcm_data(i).sigRes = 1;

            else

                vcm_data(i).sigRes = 0;

            end


            %% Recording-channel rs-fMRI time series

            recCh_ts = find( ...
                ismember(string({timeSeries.name}), string(recCh)), ...
                1);

            if ~isempty(recCh_ts)

                ts = ...
                    timeSeries(recCh_ts).ts;

                vcm_data(i).meanTimeSeries = ...
                    mean(ts, 1);

                clear ts

            end

        end


        %% ----------------------------------------------------------------
        % Save integrated data for this stimulation pair
        % -----------------------------------------------------------------

        outputDir = fullfile( ...
            derivativesPath, ...
            'vcm_data', ...
            ['sub-' subj]);

        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        saveName = fullfile( ...
            outputDir, ...
            [subj '_' stim_pair '.mat']);

        save( ...
            saveName, ...
            'vcm_data', ...
            '-v7.3');

        clear vcm_data

    end

end


%% ------------------------------------------------------------------------
% Save updated vcm_Subjects
% -------------------------------------------------------------------------

save( ...
    subjectsFile, ...
    'vcm_Subjects', ...
    '-v7.3');


fprintf('\n');
toc;