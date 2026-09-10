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
%      relative distance < -1, and electrodes labeled exactly as SOZ.
%   6. Saves the updated loc_info table.
%   7. Identifies real bipolar stimulation pairs from CCEP events files.
%   8. Keeps a pair when at least one contact passes electrode selection.
%   9. Requires the run to contain events, channels, and MEF data.
%  10. If a pair occurs in multiple usable runs, selects the run with the
%      largest number of trials; ties are resolved using the first
%      numerical run.
%  11. Returns one row per bipolar stimulation pair in vcm_Subjects.
%
% INPUTS
%   localDataPath - Path to local BIDS dataset
%   all_subjects  - Cell array of subject labels
%
% OUTPUT
%   vcm_Subjects  - Structure containing selected bipolar stimulation pairs
%                   by subject
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
            ['Electrode order differs between electrode and ' ...
             'GM/WM tables for sub-%s.'], ...
            sub_label);
    end

    if ~isequal( ...
            string(loc_info.name), ...
            string(yeoInfo.Electrode_name))

        error( ...
            ['Electrode order differs between electrode and ' ...
             'Yeo tables for sub-%s.'], ...
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
    % Reproduces historical Yeo7SubjMaxVertices behavior:
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

    hasCCEP = ...
        loc_info.CCEPs > 0;

    isNetworkROI = ...
        ismember( ...
            string(loc_info.yeo7VisualControl), ...
            ["Visual", "Control"]);

    isGrayMatterEligible = ...
        loc_info.gm_wm_relativeDistance >= -1;

    % Preserve historical exact-match behavior.
    isNotSOZ = ...
        string(loc_info.seizure_zone) ~= "SOZ";

    Ch2work = find( ...
        hasCCEP & ...
        isNetworkROI & ...
        isGrayMatterEligible & ...
        isNotSOZ);

    selectedElectrodes = ...
        string(loc_info.name(Ch2work));

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

    %% Find real bipolar stimulation pairs

    ieegDir = fullfile( ...
        localDataPath, ...
        ['sub-' sub_label], ...
        ['ses-' ses_label], ...
        'ieeg');

    eventFiles = dir(fullfile( ...
        ieegDir, ...
        '*_task-ccep_*events.tsv'));

    if isempty(eventFiles)
        warning( ...
            'No CCEP events files found for sub-%s.', ...
            sub_label);
    end

    % Store all usable pair/run combinations.
    allPairs = strings(0, 1);
    allRuns = strings(0, 1);
    allTrials = zeros(0, 1);

    allEl1 = strings(0, 1);
    allEl2 = strings(0, 1);

    allPairNetworks = strings(0, 1);

    nNetworkConflicts = 0;

    for ff = 1:numel(eventFiles)

        eventsFile = fullfile( ...
            eventFiles(ff).folder, ...
            eventFiles(ff).name);

        %% Identify run

        runNumber = extractBetween( ...
            string(eventFiles(ff).name), ...
            'ccep_run-', ...
            '_events');

        if isempty(runNumber) || ismissing(runNumber)

            warning( ...
                'Could not identify run number from %s.', ...
                eventFiles(ff).name);

            continue

        end

        %% Require files needed for BSEP preprocessing

        fileStem = [ ...
            'sub-' sub_label ...
            '_ses-' ses_label ...
            '_task-ccep_run-' char(runNumber)];

        channelsFile = fullfile( ...
            ieegDir, ...
            [fileStem '_channels.tsv']);

        mefFile = fullfile( ...
            ieegDir, ...
            [fileStem '_ieeg.mefd']);

        % A run is usable only if channel metadata and MEF data exist.
        if ~isfile(channelsFile) || exist(mefFile, 'dir') ~= 7

            fprintf( ...
                ['Skipping run-%s for sub-%s: ' ...
                 'preprocessing data incomplete.\n'], ...
                runNumber, ...
                sub_label);

            continue

        end

        %% Load events

        ev = readtable( ...
            eventsFile, ...
            'FileType', 'text', ...
            'Delimiter', '\t');

        if ~ismember( ...
                'electrical_stimulation_site', ...
                ev.Properties.VariableNames)

            warning( ...
                ['Skipping events file without ' ...
                 'electrical_stimulation_site: %s'], ...
                eventFiles(ff).name);

            continue

        end

        stimSites = ...
            string(ev.electrical_stimulation_site);

        uniquePairs = ...
            unique(stimSites, 'stable');

        %% Evaluate stimulation pairs

        for pp = 1:numel(uniquePairs)

            stimPair = uniquePairs(pp);

            if ismissing(stimPair) || strlength(stimPair) == 0
                continue
            end

            pairParts = split(stimPair, '-');

            if numel(pairParts) ~= 2
                continue
            end

            el1 = strtrim(pairParts(1));
            el2 = strtrim(pairParts(2));

            el1Selected = ...
                ismember(el1, selectedElectrodes);

            el2Selected = ...
                ismember(el2, selectedElectrodes);

            % Keep pair if at least one contact passed selection.
            if ~(el1Selected || el2Selected)
                continue
            end

            %% Count trials for this exact pair in this run

            pairIdx = ...
                stimSites == stimPair;

            nTrials = ...
                sum(pairIdx);

            if nTrials <= 1
                continue
            end

            %% Get network assignments

            idx1 = find( ...
                string(loc_info.name) == el1, ...
                1);

            idx2 = find( ...
                string(loc_info.name) == el2, ...
                1);

            net1 = "";
            net2 = "";

            if ~isempty(idx1)
                net1 = string( ...
                    loc_info.yeo7VisualControl{idx1});
            end

            if ~isempty(idx2)
                net2 = string( ...
                    loc_info.yeo7VisualControl{idx2});
            end

            %% Assign network to stimulation pair

            pairNetwork = "";

            if el1Selected && ~el2Selected

                pairNetwork = net1;

            elseif ~el1Selected && el2Selected

                pairNetwork = net2;

            elseif el1Selected && el2Selected

                if net1 == net2

                    pairNetwork = net1;

                else

                    nNetworkConflicts = ...
                        nNetworkConflicts + 1;

                    warning( ...
                        ['Skipping cross-network stimulation pair %s ' ...
                         'in sub-%s (%s vs %s).'], ...
                        stimPair, ...
                        sub_label, ...
                        net1, ...
                        net2);

                    continue

                end

            end

            %% Store pair/run combination

            allPairs(end+1, 1) = ...
                stimPair;

            allRuns(end+1, 1) = ...
                runNumber;

            allTrials(end+1, 1) = ...
                nTrials;

            allEl1(end+1, 1) = ...
                el1;

            allEl2(end+1, 1) = ...
                el2;

            allPairNetworks(end+1, 1) = ...
                pairNetwork;

        end

    end

    %% Select one run per bipolar stimulation pair
    %
    % If the same pair occurs in multiple usable runs:
    %
    %   1. Prefer the run with the largest number of trials.
    %   2. If runs have the same number of trials, select the first
    %      numerical run.

    if isempty(allPairs)

        stimEl1 = strings(0, 1);
        stimEl2 = strings(0, 1);
        run = cell(0, 1);
        networkYeo = strings(0, 1);

    else

        uniqueStimPairs = ...
            unique(allPairs, 'stable');

        keepIdx = ...
            zeros(numel(uniqueStimPairs), 1);

        for ii = 1:numel(uniqueStimPairs)

            idx = find( ...
                allPairs == uniqueStimPairs(ii));

            %% Criterion 1: largest number of trials

            pairTrials = ...
                allTrials(idx);

            maxTrials = ...
                max(pairTrials);

            bestIdx = ...
                idx(pairTrials == maxTrials);

            %% Criterion 2: first numerical run

            if numel(bestIdx) > 1

                runNumbers = ...
                    str2double(allRuns(bestIdx));

                if all(~isnan(runNumbers))

                    [~, jj] = ...
                        min(runNumbers);

                else

                    % Fallback for unexpected non-numeric run labels.
                    [~, order] = ...
                        sort(allRuns(bestIdx));

                    jj = ...
                        order(1);

                end

                keepIdx(ii) = ...
                    bestIdx(jj);

            else

                keepIdx(ii) = ...
                    bestIdx;

            end

        end

        %% Final pair-level variables

        stimEl1 = ...
            allEl1(keepIdx);

        stimEl2 = ...
            allEl2(keepIdx);

        selectedRuns = ...
            allRuns(keepIdx);

        networkYeo = ...
            allPairNetworks(keepIdx);

        % Preserve historical vcm_Subjects representation:
        % each table cell contains one run string.
        run = cell(numel(selectedRuns), 1);

        for ii = 1:numel(selectedRuns)
            run{ii} = selectedRuns(ii);
        end

    end

    %% Create subject stimulation-pair table

    subjField = ...
        matlab.lang.makeValidName(sub_label);

    vcm_Subjects.(subjField) = table( ...
        stimEl1, ...
        stimEl2, ...
        run, ...
        networkYeo);

    %% Print summary

    fprintf( ...
        'sub-%s: %d electrodes passed selection criteria.\n', ...
        sub_label, ...
        numel(Ch2work));

    fprintf( ...
        'sub-%s: %d bipolar stimulation pairs selected.\n', ...
        sub_label, ...
        height(vcm_Subjects.(subjField)));

    if ~isempty(networkYeo)

        [networkNames, ~, groupIdx] = ...
            unique(networkYeo, 'stable');

        networkCounts = ...
            accumarray(groupIdx, 1);

        for ii = 1:numel(networkNames)

            fprintf( ...
                '  %s pairs: %d\n', ...
                networkNames(ii), ...
                networkCounts(ii));

        end

    end

    if nNetworkConflicts > 0

        fprintf( ...
            '  Cross-network pairs skipped: %d\n', ...
            nNetworkConflicts);

    end

end

fprintf('\n');
toc;

end