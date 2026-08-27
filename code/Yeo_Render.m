%% plot inflated brains
function Yeo_Render(localDataPath, all_subjects, hem, yeo, targetNetwork)
% Yeo_Render
%
% Purpose:
%   Plot Yeo networks on the subject inflated surface with electrodes.
%
% Optional:
%   targetNetwork = 'all'          -> plot all Yeo networks
%   targetNetwork = 'Visual'       -> plot only Visual network
%   targetNetwork = 'Somatomotor'  -> plot only Somatomotor network
%   targetNetwork = 'Control'      -> plot only Control network
%
% Example:
%   Yeo_Render(localDataPath, {'01'}, 1, 7, 'Visual')
%   Yeo_Render(localDataPath, {'01'}, [1 2], 7, 'all')
%
% IMPORTANT:
%   First run the Yeo surf2surf step so the subject has:
%     lh.Yeo2011_7Networks_N1000.annot
%     rh.Yeo2011_7Networks_N1000.annot
%
% Inputs:
%   localDataPath = project data folder
%   all_subjects  = cell array of subjects, e.g. {'01'}
%   hem           = 1 for left, 2 for right, or [1 2]
%   yeo           = 7 or 17
%   targetNetwork = optional network name or 'all'

    if nargin < 5 || isempty(targetNetwork)
        targetNetwork = 'all';
    end

    for ss = 1:length(all_subjects)

        sub_label = all_subjects{ss};
        ses_label = 'ieeg01';

        fprintf('\nRendering subject: %s\n', sub_label);
        fprintf('Yeo: %d\n', yeo);
        fprintf('Target network: %s\n', targetNetwork);

        %% get elecmatrix

        electrodes_tsv_name = fullfile(localDataPath, ...
            'derivatives/electrodes2rsfMRI', ...
            ['sub-' sub_label], ...
            ['sub-' sub_label '_ses-' ses_label '_space-rsT1_electrodes.tsv']);

        loc_info = readtable(electrodes_tsv_name, ...
            'FileType', 'text', ...
            'Delimiter', '\t', ...
            'TreatAsEmpty', {'N/A','n/a'});

        %% load pial and inflated giftis

        gL_infl = gifti(fullfile(localDataPath, ...
            'derivatives/freesurfer', ...
            ['sub-' sub_label], ...
            'inflated.L.surf.gii'));

        gR_infl = gifti(fullfile(localDataPath, ...
            'derivatives/freesurfer', ...
            ['sub-' sub_label], ...
            'inflated.R.surf.gii'));

        gL = gifti(fullfile(localDataPath, ...
            'derivatives/freesurfer', ...
            ['sub-' sub_label], ...
            'pial.L.surf.gii'));

        gR = gifti(fullfile(localDataPath, ...
            'derivatives/freesurfer', ...
            ['sub-' sub_label], ...
            'pial.R.surf.gii'));

        %% snap electrodes to surface and then move to inflated

        xyz_inflated = ieeg_snap2inflated(loc_info, gR, gL, gR_infl, gL_infl);

        views_plot = {[65, -50], [-60, 30]};

        %% loop over hemisphere

        for hh = hem

            if hh == 1
                hemi = 'l';
                g = gL_infl;
            elseif hh == 2
                hemi = 'r';
                g = gR_infl;
            else
                error('hem should be 1, 2, or [1 2]');
            end

            %% Surface labels

            if yeo == 7

                [~, verlab, temp] = read_annotation(fullfile(localDataPath, ...
                    'derivatives', ...
                    'freesurfer', ...
                    ['sub-' sub_label], ...
                    'label', ...
                    [hemi 'h.Yeo2011_7Networks_N1000.annot']));

                % Yeo 7 labeling
                Yeo_ROI_Names = {'Visual';          % 1
                                 'Somatomotor';     % 2
                                 'DorsalAttention'; % 3
                                 'VentralAttention';% 4
                                 'Limbic';          % 5
                                 'Control';         % 6
                                 'Default'};        % 7

            elseif yeo == 17

                [~, verlab, temp] = read_annotation(fullfile(localDataPath, ...
                    'derivatives', ...
                    'freesurfer', ...
                    ['sub-' sub_label], ...
                    'label', ...
                    [hemi 'h.Yeo2011_17Networks_N1000.annot']));

                % Yeo 17 labeling
                Yeo_ROI_Names = {'VisualCentral';      % 1
                                 'VisualPeripheral';   % 2
                                 'SomatomotorA';       % 3
                                 'SomatomotorB';       % 4
                                 'DorsalAttentionA';   % 5
                                 'DorsalAttentionB';   % 6
                                 'VentralAttentionA';  % 7
                                 'VentralAttentionB';  % 8
                                 'LimbicA';            % 9
                                 'LimbicB';            % 10
                                 'ControlC';           % 11
                                 'ControlA';           % 12
                                 'ControlB';           % 13
                                 'TemporalParietal';   % 14
                                 'DefaultC';           % 15
                                 'DefaultA';           % 16
                                 'DefaultB'};          % 17

            else
                error('Yeo should be 7 or 17');
            end

            %% Yeo colors

            temp.table(1:size(temp.table,1),6) = 0:1:(size(temp.table,1)-1);

            cmap = temp.table(2:end, 1:3) ./ 255;
            cmap = cmap.^0.4;

            for verS = 1:length(temp.table)
                idxCh = find(ismember(verlab, temp.table(verS,5)));
                verlab(idxCh) = temp.table(verS,6);
            end

            vert_label = verlab;

           %% Optional: plot only one network

            if ~strcmpi(targetNetwork, 'all')
            
                networkIdx = find(strcmpi(Yeo_ROI_Names, targetNetwork));
            
                if isempty(networkIdx)
                    error('Network "%s" was not found. Check spelling.', targetNetwork);
                end
            
                % Save original color for selected network
                targetColor = cmap(networkIdx, :);
            
                % Keep only the selected network.
                % Everything else becomes NaN, so it should not be colored.
                selected_label = NaN(size(vert_label));
                selected_label(vert_label == networkIdx) = 1;
            
                vert_label = selected_label;
            
                % Only one color: the selected network color
                cmap = targetColor;
            
                Yeo_ROI_Names = {Yeo_ROI_Names{networkIdx}};
            
            end

            %% sulcal labels individuals

            sulcal_labels = read_curv(fullfile(localDataPath, ...
                'derivatives/freesurfer', ...
                ['sub-' sub_label], ...
                'surf', ...
                [hemi 'h.sulc']));

            electrodes_thisHemi = find(ismember(loc_info.hemisphere, upper(hemi)));

            %% make a plot with electrode dots

            for vv = 1:length(views_plot)

                v_d = [views_plot{vv}(1), views_plot{vv}(2)];

                % get the inflated coordinates
                els = xyz_inflated;

                % with labels - popout
                a_offset = .1 * max(abs(els(:,1))) * ...
                    [cosd(v_d(1)-90)*cosd(v_d(2)), ...
                     sind(v_d(1)-90)*cosd(v_d(2)), ...
                     sind(v_d(2))];

                els_pop = els + repmat(a_offset, size(els,1), 1);

                figure

                tH = ieeg_RenderGiftiLabels( ...
                    g, ...
                    vert_label, ...
                    cmap, ...
                    Yeo_ROI_Names, ...
                    sulcal_labels); %#ok<NASGU>

                % ieeg_elAdd( ...
                %     els_pop(electrodes_thisHemi,:), ...
                %     [.1 .1 .1], ...
                %     15); % add electrode positions

                % Optional: add electrode names
                % ieeg_label( ...
                %     els_pop(electrodes_thisHemi,:), ...
                %     15, ...
                %     10, ...
                %     loc_info.name(electrodes_thisHemi));

                ieeg_viewLight(v_d(1), v_d(2)); % change viewing angle

                title([sub_label ' - ' upper(hemi) 'H - Yeo' num2str(yeo) ...
                    ' - ' targetNetwork], ...
                    'Interpreter', 'none');

            end

        end

    end

end