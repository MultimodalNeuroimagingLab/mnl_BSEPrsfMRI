function vcm_Subjects = ssf_select_subject_electrodes(localDataPath, all_subjects)
% ssf_select_subject_electrodes
%
% Select stimulation electrodes with available CCEP/BSEP data located in
% the Yeo Visual or Control networks and satisfying anatomical criteria.
%
% For each subject, this function:
%   1. Loads electrode coordinates in rs-fMRI T1 space.
%   2. Identifies electrodes with available CCEP stimulation trials.
%   3. Assigns each electrode to Visual or Control based on the number of
%      nearby Yeo-7 surface vertices.
%   4. Adds gray/white matter relative distance.
%   5. Excludes electrodes outside the selected networks, electrodes with
%      relative distance < -1, and electrodes labeled as SOZ.
%   6. Saves the updated loc_info table.
%   7. Returns the selected electrodes in vcm_Subjects.
%
% INPUTS
%   localDataPath - Path to local BIDS dataset
%   all_subjects  - Cell array of subject labels
%
% OUTPUT
%   vcm_Subjects  - Structure containing selected electrodes by subject
%
% Author: Maria Guadalupe Yanez Ramos
% Developed with scientific and technical guidance from Dora Hermes
% and the Multimodal Neuroimaging Lab (MNL) team.
% May 2026

tic;

ses_label = 'ieeg01';

vcm_Subjects = struct();

for ss = 1:numel(all_subjects)

    sub_label = all_subjects{ss};

    fprintf('\nProcessing sub-%s...\n', sub_label);

    %% Load transformed electrode information

    electrodesFile = fullfile( ...
        localDataPath, ...
        'derivatives', ...
        'electrodes2rsfMRI', ...
        ['sub-' sub_label], ...
        ['sub-' sub_label '_ses-' ses_label ...
         '_space-rsT1_electrodes.tsv']);

    loc_info = readtable( ...
        electrodesFile, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'TreatAsEmpty', {'N/A', 'n/a'});

    nElectrodes = height(loc_info);

    %% Load gray/white matter distance

    gmFile = fullfile( ...
        localDataPath, ...
        'derivatives', ...
        'gm_wmDistance', ...
        ['sub-' sub_label], ...
        ['sub-' sub_label '_ses-' ses_label '_gm_wm.tsv']);

    gm_wm_info = readtable( ...
        gmFile, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'TreatAsEmpty', {'N/A', 'n/a'});

    %% Load Yeo-7 electrode information

    yeoFile = fullfile( ...
        localDataPath, ...
        'derivatives', ...
        'Yeo7Electrodes', ...
        ['sub-' sub_label '_TableYeoElectrodes.xlsx']);

    yeoInfo = readtable(yeoFile);

    %% Check electrode correspondence

    if height(gm_wm_info) ~= nElectrodes
        error( ...
            'GM/WM electrode count does not match for sub-%s.', ...
            sub_label);
    end

    if height(yeoInfo) ~= nElectrodes
        error( ...
            'Yeo electrode count does not match for sub-%s.', ...
            sub_label);
    end

    if ~isequal( ...
            string(loc_info.name), ...
            string(gm_wm_info.ElectrodeName))

        error( ...
            'Electrode order differs between electrode and GM/WM tables for sub-%s.', ...
            sub_label);
    end

    if ~isequal( ...
            string(loc_info.name), ...
            string(yeoInfo.Electrode_name))

        error( ...
            'Electrode order differs between electrode and Yeo tables for sub-%s.', ...
            sub_label);
    end

    %% Initialize variables added to loc_info

    loc_info.CCEPs = zeros(nElectrodes, 1);
    loc_info.runCCEPS = cell(nElectrodes, 1);
    loc_info.yeo7VisualControl = cell(nElectrodes, 1);

    loc_info.gm_wm_relativeDistance = ...
        gm_wm_info.RelativeDistance;

    %% Determine CCEP availability

    for eli = 1:nElectrodes

        Ch = char(string(loc_info.name(eli)));
        [ccepExist, runccepExist] = ...
            ssf_find_ccep_thisCh( ...
                Ch, ...
                sub_label, ...
                localDataPath);

        % Historical criterion:
        % require more than one stimulation trial.
        if ccepExist > 1

            loc_info.CCEPs(eli) = 1;

            loc_info.runCCEPS{eli} = ...
                extractBetween( ...
                    runccepExist, ...
                    'ccep_run-', ...
                    '_event');

        else

            loc_info.CCEPs(eli) = 0;
            loc_info.runCCEPS{eli} = [];

        end

    end

    %% Assign Visual or Control network
    %
    % Reproduces the historical Yeo7SubjMaxVertices behavior:
    %
    %   Control > Visual     -> Control
    %   Visual = Control = 0 -> no network
    %   Visual >= Control    -> Visual

    for eli = 1:nElectrodes

        visualVertices = yeoInfo.Visual(eli);
        controlVertices = yeoInfo.Control(eli);

        if visualVertices < controlVertices

            loc_info.yeo7VisualControl{eli} = 'Control';

        elseif visualVertices == 0 && controlVertices == 0

            loc_info.yeo7VisualControl{eli} = '';

        else

            loc_info.yeo7VisualControl{eli} = 'Visual';

        end

    end

    %% Select electrodes

    hasCCEP = loc_info.CCEPs > 0;

    isNetworkROI = ...
        ismember(string(loc_info.yeo7VisualControl), ...
        ["Visual", "Control"]);

    isGrayMatterEligible = ...
        loc_info.gm_wm_relativeDistance >= -1;

    isNotSOZ = ...
        string(loc_info.seizure_zone) ~= "SOZ";

    Ch2work = find( ...
        hasCCEP & ...
        isNetworkROI & ...
        isGrayMatterEligible & ...
        isNotSOZ);

    %% Save updated loc_info

    locInfoOutputDir = fullfile( ...
        localDataPath, ...
        'derivatives', ...
        'loc_info', ...
        ['sub-' sub_label]);

    if ~exist(locInfoOutputDir, 'dir')
        mkdir(locInfoOutputDir);
    end

    locInfoFile = fullfile( ...
        locInfoOutputDir, ...
        'loc_info.mat');

    save(locInfoFile, 'loc_info');

    %% Create subject selection table

    subjField = matlab.lang.makeValidName(sub_label);

    stimEl1 = string(loc_info.name(Ch2work));

    % Historical behavior preserved temporarily.
    % stimEl2 must later be replaced by the actual second electrode
    % of the bipolar stimulation pair before BSEP preprocessing.
    stimEl2 = string(loc_info.name(Ch2work));

    run = loc_info.runCCEPS(Ch2work);

    networkYeo = ...
        string(loc_info.yeo7VisualControl(Ch2work));

    vcm_Subjects.(subjField) = table( ...
        stimEl1, ...
        stimEl2, ...
        run, ...
        networkYeo);

    %% Print summary

    fprintf( ...
        'sub-%s: %d electrodes selected.\n', ...
        sub_label, ...
        numel(Ch2work));

    if ~isempty(Ch2work)

        selectedNetworks = ...
            string(loc_info.yeo7VisualControl(Ch2work));

        [networkNames, ~, groupIdx] = ...
            unique(selectedNetworks, 'stable');

        networkCounts = ...
            accumarray(groupIdx, 1);

        for ii = 1:numel(networkNames)

            fprintf( ...
                '  %s: %d\n', ...
                networkNames(ii), ...
                networkCounts(ii));

        end

    end

end

fprintf('\n');
toc;

end