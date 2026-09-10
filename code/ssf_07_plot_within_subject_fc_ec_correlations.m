%% ssf_07_plot_within_subject_fc_ec_correlations
%
% Plot within-subject relationships between resting-state functional
% connectivity (FC Fisher z) and effective connectivity (EC CoD).
%
% Columns correspond to stimulation network.
%
% Top row:
%   recording channels within the stimulation network.
%
% Bottom row:
%   all valid recording channels.
%
% Individual subjects are shown with separate colors. Regression lines
% are estimated independently within each subject.
%
% INPUT
%   data/allRows.mat
%
% OUTPUT
%   derivatives/figures/fc_ec_correlations/
%
% Author: Maria Guadalupe Yanez Ramos
% Developed with scientific and technical guidance from Dora Hermes
% and the Multimodal Neuroimaging Lab (MNL) team.
% May 2026


clearvars;
close all;
clc;


%% ------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------

scriptPath = which( ...
    'ssf_07_plot_within_subject_fc_ec_correlations');

if isempty(scriptPath)
    error('Could not determine the Stage 07 script location.');
end

codeDir = fileparts(scriptPath);
projectDir = fileparts(codeDir);

dataPath = fullfile(projectDir, 'data');

figuresPath = fullfile( ...
    dataPath, ...
    'derivatives', ...
    'figures', ...
    'fc_ec_correlations');

if ~exist(figuresPath, 'dir')
    mkdir(figuresPath);
end


%% ------------------------------------------------------------------------
% Load Stage 06 table
% -------------------------------------------------------------------------

S = load( ...
    fullfile(dataPath, 'allRows.mat'), ...
    'allRows');

allRows = S.allRows;


%% ------------------------------------------------------------------------
% Keep observations with valid FC and EC
% -------------------------------------------------------------------------

T = allRows( ...
    allRows.Add2Analyses == 1 & ...
    isfinite(allRows.zFisher) & ...
    isfinite(allRows.mean_cod_noOutliers), ...
    :);


% Convert variables used repeatedly below.

T.subject = string(T.subject);
T.stim_network = string(T.stim_network);


%% ------------------------------------------------------------------------
% Determine stimulation networks represented in current data
% -------------------------------------------------------------------------

preferredOrder = ...
    ["Visual", "Somatomotor", "Control"];

presentNetworks = ...
    unique(T.stim_network, 'stable');

networks = ...
    preferredOrder(ismember(preferredOrder, presentNetworks));


if isempty(networks)
    error('No supported stimulation networks were found.');
end


%% ------------------------------------------------------------------------
% Subjects and colors
% -------------------------------------------------------------------------

subjects = ...
    unique(T.subject, 'stable');

subjectColors = ...
    lines(numel(subjects));


%% ------------------------------------------------------------------------
% Figure
% -------------------------------------------------------------------------

fig = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [100 100 420*numel(networks) 720]);


tl = tiledlayout( ...
    2, ...
    numel(networks), ...
    'TileSpacing', ...
    'compact', ...
    'Padding', ...
    'compact');


%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

for row = 1:2

    for n = 1:numel(networks)

        net = networks(n);

        ax = nexttile;
        hold(ax, 'on');


        %% ---------------------------------------------------------------
        % Select panel data
        % ----------------------------------------------------------------

        if row == 1

            % Recording channel belongs to the stimulated network.

            Tpanel = T( ...
                T.stim_network == net & ...
                T.insideStimNetwork == 1, ...
                :);

        else

            % All valid recording channels.

            Tpanel = T( ...
                T.stim_network == net, ...
                :);

        end


        %% ---------------------------------------------------------------
        % Plot each subject independently
        % ----------------------------------------------------------------

        for s = 1:numel(subjects)

            subj = subjects(s);

            Tsub = Tpanel( ...
                Tpanel.subject == subj, ...
                :);

            if isempty(Tsub)
                continue
            end


            x = double(Tsub.zFisher(:));
            y = double(Tsub.mean_cod_noOutliers(:));

            good = ...
                isfinite(x) & ...
                isfinite(y);

            x = x(good);
            y = y(good);


            if isempty(x)
                continue
            end


            %% Points

            scatter( ...
                ax, ...
                x, ...
                y, ...
                14, ...
                subjectColors(s,:), ...
                'filled', ...
                'MarkerFaceAlpha', ...
                0.25, ...
                'MarkerEdgeAlpha', ...
                0.25, ...
                'DisplayName', ...
                char(subj));


            %% Within-subject regression line

            if numel(x) >= 3 && ...
                    std(x) > 0 && ...
                    std(y) > 0

                lm = ...
                    fitlm(x, y);

                xFit = ...
                    linspace( ...
                        min(x), ...
                        max(x), ...
                        100)';

                yFit = ...
                    predict(lm, xFit);

                plot( ...
                    ax, ...
                    xFit, ...
                    yFit, ...
                    '-', ...
                    'Color', ...
                    subjectColors(s,:), ...
                    'LineWidth', ...
                    2.2, ...
                    'HandleVisibility', ...
                    'off');

            end

        end


        %% ---------------------------------------------------------------
        % Panel formatting
        % ----------------------------------------------------------------

        if row == 1

            title( ...
                ax, ...
                sprintf('%s: within-network', net));
            xlim([0 2.5]);
            ylim([-0.5 1]);

        else

            title( ...
                ax, ...
                sprintf('%s: all recording sites', net));
            xlim([0 2.5]);
            ylim([-0.5 1]);

        end


        if row == 2
            xlabel(ax, 'FC (Fisher z)');
            xlim([0 2.5]);
            ylim([-0.5 1]);
        end


        if n == 1
            ylabel(ax, 'EC (median CoD)');
            xlim([0 2.5]);
            ylim([-0.5 1]);
        end


        box(ax, 'off');

        set( ...
            ax, ...
            'FontSize', ...
            11);

    end

end


title( ...
    tl, ...
    'Correlations between FC and EC within subjects', ...
    'FontSize', ...
    18, ...
    'FontWeight', ...
    'bold');


%% ------------------------------------------------------------------------
% Legend
% -------------------------------------------------------------------------

lgd = legend( ...
    tl.Children(end), ...
    'Location', ...
    'best');

lgd.Title.String = ...
    'Subject';


%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

outputPng = fullfile( ...
    figuresPath, ...
    'ssf_within_subject_fc_ec_correlations.png');

outputFig = fullfile( ...
    figuresPath, ...
    'ssf_within_subject_fc_ec_correlations.fig');


exportgraphics( ...
    fig, ...
    outputPng, ...
    'Resolution', ...
    300);

savefig( ...
    fig, ...
    outputFig);


fprintf('\nValid rows: %d\n', height(T));

fprintf('\nRows by stimulation network:\n');
disp(groupsummary(T, 'stim_network'));

fprintf('\nWithin-network rows:\n');

Twithin = T( ...
    T.insideStimNetwork == 1, ...
    :);

disp(groupsummary( ...
    Twithin, ...
    'stim_network'));

fprintf('\nSaved figure:\n%s\n', outputPng);