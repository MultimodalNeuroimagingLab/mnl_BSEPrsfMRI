%% ssf_06_compute_fc_analysis_table
%
% Compute resting-state functional connectivity between each stimulation
% seed and recording-channel ROI, and assemble the analysis table used for
% FC-BSEP correlation analyses.
%
% INPUTS
%   data/vcm_Subjects.mat
%   derivatives/vcm_data/sub-<label>/*.mat
%
% OUTPUTS
%   updated derivatives/vcm_data/sub-<label>/*.mat
%   data/allRows.mat
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

scriptPath = which('ssf_06_compute_fc_analysis_table');

if isempty(scriptPath)
    error( ...
        ['Could not determine the location of ' ...
         'ssf_06_compute_fc_analysis_table.m.']);
end

codeDir = fileparts(scriptPath);
projectDir = fileparts(codeDir);

localDataPath = fullfile(projectDir, 'data');
derivativesPath = fullfile(localDataPath, 'derivatives');


%% ------------------------------------------------------------------------
% Load subject/stimulation-pair information
% -------------------------------------------------------------------------

subjectsFile = fullfile( ...
    localDataPath, ...
    'vcm_Subjects.mat');

if ~isfile(subjectsFile)
    error('vcm_Subjects.mat was not found.');
end

tmp = load(subjectsFile, 'vcm_Subjects');
vcm_Subjects = tmp.vcm_Subjects;

all_subjects = fieldnames(vcm_Subjects)';


%% ------------------------------------------------------------------------
% Analysis settings
% -------------------------------------------------------------------------

% Historical recording-channel inclusion threshold used in vcm06.
gmWmThreshold = -1.4;


%% ------------------------------------------------------------------------
% Initialize long-format analysis structure
% -------------------------------------------------------------------------

rowN = 0;

rows = struct( ...
    'subject', {}, ...
    'stim_pair', {}, ...
    'stim_network', {}, ...
    'rec_channel', {}, ...
    'rec_network', {}, ...
    'insideStimNetwork', {}, ...
    'gm_wm_relativeDistance', {}, ...
    'Add2Analyses', {}, ...
    'rPearson', {}, ...
    'pPearson', {}, ...
    'zFisher', {}, ...
    'median_cod', {}, ...
    'mean_cod_noOutliers', {}, ...
    'sigRes', {});


%% ------------------------------------------------------------------------
% Process subjects
% -------------------------------------------------------------------------

for iSubj = 1:numel(all_subjects)

    subj = all_subjects{iSubj};

    subjectInfo = ...
        vcm_Subjects.(subj);

    nStimPairs = ...
        height(subjectInfo);

    fprintf( ...
        '\nProcessing sub-%s (%d of %d)\n', ...
        subj, ...
        iSubj, ...
        numel(all_subjects));


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

        stimNetwork = ...
            string(subjectInfo.networkYeo(j));

        fprintf( ...
            'Stim pair %d of %d: %s (%s)\n', ...
            j, ...
            nStimPairs, ...
            stim_pair, ...
            stimNetwork);


        %% ----------------------------------------------------------------
        % Load integrated BSEP/rs-fMRI data from Stage 05
        % -----------------------------------------------------------------

        vcmFile = fullfile( ...
            derivativesPath, ...
            'vcm_data', ...
            ['sub-' subj], ...
            [subj '_' stim_pair '.mat']);

        if ~isfile(vcmFile)
            error( ...
                'vcm_data file not found for sub-%s, pair %s.', ...
                subj, ...
                stim_pair);
        end

        tmpVcm = load(vcmFile, 'vcm_data');
        vcm_data = tmpVcm.vcm_data;


        %% ----------------------------------------------------------------
        % Stimulation seed
        % -----------------------------------------------------------------

        seed = ...
            subjectInfo.meanTimeSeries{j};

        if isempty(seed)
            error( ...
                'Empty stimulation seed for sub-%s, pair %s.', ...
                subj, ...
                stim_pair);
        end

        seed = ...
            double(seed(:));


        %% ----------------------------------------------------------------
        % Recording channels
        % -----------------------------------------------------------------

        for i = 1:numel(vcm_data)

            recCh = ...
                char(string(vcm_data(i).use_channels_name));


            %% ------------------------------------------------------------
            % Existing metadata
            % -------------------------------------------------------------

            recNetwork = "";

            if ~isempty(vcm_data(i).yeo7VisualControl)
                recNetwork = ...
                    string(vcm_data(i).yeo7VisualControl);
            end


            gmWmDistance = NaN;

            if ~isempty(vcm_data(i).gm_wm_relativeDistance)
                gmWmDistance = ...
                    double(vcm_data(i).gm_wm_relativeDistance);
            end


            medianCod = NaN;

            if ~isempty(vcm_data(i).median_cod)
                medianCod = ...
                    double(vcm_data(i).median_cod);
            end


            meanCod = NaN;

            if ~isempty(vcm_data(i).mean_cod_noOutliers)
                meanCod = ...
                    double(vcm_data(i).mean_cod_noOutliers);
            end


            sigRes = NaN;

            if ~isempty(vcm_data(i).sigRes)
                sigRes = ...
                    double(vcm_data(i).sigRes);
            end


            %% ------------------------------------------------------------
            % Initialize FC values
            % -------------------------------------------------------------

            addToAnalysis = 0;

            rPearson = NaN;
            pPearson = NaN;
            zFisher = NaN;


            %% ------------------------------------------------------------
            % Compute rs-fMRI FC
            % -------------------------------------------------------------

            hasRsFmri = ...
                ~isempty(vcm_data(i).meanTimeSeries);

            validDistance = ...
                isfinite(gmWmDistance) && ...
                gmWmDistance >= gmWmThreshold;


            if hasRsFmri && validDistance

                roi = ...
                    double(vcm_data(i).meanTimeSeries(:));

                if numel(roi) < numel(seed)

                    error( ...
                        ['Recording-channel time series is shorter than ' ...
                         'the stimulation seed for sub-%s, pair %s, ' ...
                         'channel %s.'], ...
                        subj, ...
                        stim_pair, ...
                        recCh);

                end


                % The stimulation seed had the first 10 volumes removed
                % during Stage 05. Align the recording-channel signal to
                % the same final volumes.

                volStart = ...
                    numel(roi) - numel(seed) + 1;

                roiAligned = ...
                    roi(volStart:end);


                % Pearson correlation is undefined for a constant signal.
                if std(seed) > 0 && std(roiAligned) > 0

                    [rPearson, pPearson] = ...
                        corr( ...
                            seed, ...
                            roiAligned, ...
                            'Type', ...
                            'Pearson');

                    zFisher = ...
                        atanh(rPearson);


                    % Include only valid finite FC estimates.
                    if isfinite(zFisher)
                        addToAnalysis = 1;
                    end

                end

            end

            %% ------------------------------------------------------------
            % Save FC results back into vcm_data
            % -------------------------------------------------------------

            vcm_data(i).Add2Analyses = ...
                addToAnalysis;

            vcm_data(i).rPearson = ...
                rPearson;

            vcm_data(i).pPearson = ...
                pPearson;

            vcm_data(i).zFisher = ...
                zFisher;


            % Same-network membership is assigned after all rows are assembled,
            % using independent Yeo network memberships.
            insideStimNetwork = NaN;


            %% ------------------------------------------------------------
            % Add one row to analysis table
            % -------------------------------------------------------------

            rowN = rowN + 1;

            rows(rowN).subject = ...
                string(subj);

            rows(rowN).stim_pair = ...
                string(stim_pair);

            rows(rowN).stim_network = ...
                stimNetwork;

            rows(rowN).rec_channel = ...
                string(recCh);

            rows(rowN).rec_network = ...
                recNetwork;

            rows(rowN).insideStimNetwork = ...
                insideStimNetwork;

            rows(rowN).gm_wm_relativeDistance = ...
                gmWmDistance;

            rows(rowN).Add2Analyses = ...
                addToAnalysis;

            rows(rowN).rPearson = ...
                rPearson;

            rows(rowN).pPearson = ...
                pPearson;

            rows(rowN).zFisher = ...
                zFisher;

            rows(rowN).median_cod = ...
                medianCod;

            rows(rowN).mean_cod_noOutliers = ...
                meanCod;

            rows(rowN).sigRes = ...
                sigRes;

        end


        %% ----------------------------------------------------------------
        % Save updated vcm_data
        % -----------------------------------------------------------------

        save( ...
            vcmFile, ...
            'vcm_data', ...
            '-v7.3');

        clear vcm_data

    end

end


%% ------------------------------------------------------------------------
% Convert to long-format table
% -------------------------------------------------------------------------

allRows = ...
    struct2table(rows);

%% ------------------------------------------------------------------------
% Independent Yeo network membership for recording contacts
% -------------------------------------------------------------------------

allRows.Somatomotor = nan(height(allRows), 1);
allRows.Visual      = nan(height(allRows), 1);
allRows.Control     = nan(height(allRows), 1);

subjects = unique(string(allRows.subject), 'stable');

for s = 1:numel(subjects)

    subj = subjects(s);

    yeoFile = fullfile( ...
        localDataPath, ...
        'derivatives', ...
        'Yeo7Electrodes', ...
        ['sub-' char(subj) '_TableYeoElectrodes.xlsx']);

    if ~isfile(yeoFile)
        error( ...
            'Yeo electrode table not found for sub-%s:\n%s', ...
            subj, yeoFile);
    end

    Y = readtable(yeoFile);

    rows = find(string(allRows.subject) == subj);

    for k = rows(:)'

        recCh = string(allRows.rec_channel(k));

        ii = find( ...
            string(Y.Electrode_name) == recCh, ...
            1);

        if isempty(ii)
            continue
        end

        allRows.Somatomotor(k) = ...
            double(Y.Somatomotor(ii) > 0);

        allRows.Visual(k) = ...
            double(Y.Visual(ii) > 0);

        allRows.Control(k) = ...
            double(Y.Control(ii) > 0);
    end
end


%% ------------------------------------------------------------------------
% Recording contact inside/outside stimulated network
% -------------------------------------------------------------------------

allRows.insideStimNetwork(:) = NaN;

idx = string(allRows.stim_network) == "Somatomotor";
allRows.insideStimNetwork(idx) = ...
    allRows.Somatomotor(idx);

idx = string(allRows.stim_network) == "Visual";
allRows.insideStimNetwork(idx) = ...
    allRows.Visual(idx);

idx = string(allRows.stim_network) == "Control";
allRows.insideStimNetwork(idx) = ...
    allRows.Control(idx);

%% ------------------------------------------------------------------------
% Basic validation summary
% -------------------------------------------------------------------------

fprintf('\nAnalysis table summary\n');
fprintf('Total rows:              %d\n', height(allRows));
fprintf('Rows Add2Analyses = 1:   %d\n', sum(allRows.Add2Analyses == 1));
fprintf('Finite FC Fisher z:       %d\n', sum(isfinite(allRows.zFisher)));
fprintf('Finite median CoD:        %d\n', sum(isfinite(allRows.median_cod)));
fprintf('Finite mean CoD:          %d\n', sum(isfinite(allRows.mean_cod_noOutliers)));

fprintf('\nRows by stimulation network:\n');
disp(groupsummary(allRows, 'stim_network'));

validRows = ...
    allRows.Add2Analyses == 1 & ...
    isfinite(allRows.zFisher) & ...
    isfinite(allRows.median_cod) & ...
    isfinite(allRows.insideStimNetwork);

fprintf('\nValid rows by stimulation and recording relation:\n');

disp(groupsummary( ...
    allRows(validRows, :), ...
    {'stim_network', 'insideStimNetwork'}));


%% ------------------------------------------------------------------------
% Save analysis table
% -------------------------------------------------------------------------

allRowsFile = fullfile( ...
    localDataPath, ...
    'allRows.mat');

save( ...
    allRowsFile, ...
    'allRows', ...
    '-v7.3');


fprintf('\nSaved:\n%s\n', allRowsFile);
toc;