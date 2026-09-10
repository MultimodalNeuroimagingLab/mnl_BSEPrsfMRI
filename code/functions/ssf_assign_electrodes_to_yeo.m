function ssf_assign_electrodes_to_yeo(localDataPath, mmDistance, all_subjects, yeo)
% ssf_assign_electrodes_to_yeo Assign Yeo network labels to electrodes.
%
% IMPORTANT: First run mri_surf2surf from fsaverage to individual space.
%
% Inputs:
%   localDataPath - Path to dataset
%   mmDistance    - Distance around electrode in mm
%   all_subjects  - Cell array of subject labels
%   yeo           - Yeo parcellation (7 or 17)
%
% Author: Maria Guadalupe Yanez Ramos
% Developed with scientific and technical guidance from Dora Hermes
% and the Multimodal Neuroimaging Lab (MNL) team.
% May 2026

    ses_label = 'ieeg01';

    for ss = 1:length(all_subjects)

        sub_label = all_subjects{ss};

        %% Load electrode coordinates

        electrodes_tsv_name = fullfile( ...
            localDataPath, ...
            'derivatives', ...
            'electrodes2rsfMRI', ...
            ['sub-' sub_label], ...
            ['sub-' sub_label '_ses-' ses_label ...
             '_space-rsT1_electrodes.tsv']);

        loc_info = readtable( ...
            electrodes_tsv_name, ...
            'FileType', 'text', ...
            'Delimiter', '\t', ...
            'TreatAsEmpty', {'N/A','n/a'});

        xyz = [loc_info.x loc_info.y loc_info.z];

        elect_num = size(xyz, 1);

        njj = 0;

        %% Load Yeo annotation

        if yeo == 7

            [~, verlabR, temp] = read_annotation(fullfile( ...
                localDataPath, ...
                'derivatives', ...
                'freesurfer', ...
                ['sub-' sub_label], ...
                'label', ...
                'rh.Yeo2011_7Networks_N1000.annot'));

            Yeo_ROI_Names = { ...
                'Visual';
                'Somatomotor';
                'DorsalAttention';
                'VentralAttention';
                'Limbic';
                'Control';
                'Default'};

            [~, verlabL, temp] = read_annotation(fullfile( ...
                localDataPath, ...
                'derivatives', ...
                'freesurfer', ...
                ['sub-' sub_label], ...
                'label', ...
                'lh.Yeo2011_7Networks_N1000.annot'));

        elseif yeo == 17

            [~, verlabR, temp] = read_annotation(fullfile( ...
                localDataPath, ...
                'derivatives', ...
                'freesurfer', ...
                ['sub-' sub_label], ...
                'label', ...
                'rh.Yeo2011_17Networks_N1000.annot'));

            Yeo_ROI_Names = { ...
                'VisualCentral';
                'VisualPeripheral';
                'SomatomotorA';
                'SomatomotorB';
                'DorsalAttentionA';
                'DorsalAttentionB';
                'VentralAttentionA';
                'VentralAttentionB';
                'LimbicA';
                'LimbicB';
                'ControlC';
                'ControlA';
                'ControlB';
                'Temporal Parietal';
                'DefaultC';
                'DefaultA';
                'DefaultB'};

            [~, verlabL, temp] = read_annotation(fullfile( ...
                localDataPath, ...
                'derivatives', ...
                'freesurfer', ...
                ['sub-' sub_label], ...
                'label', ...
                'lh.Yeo2011_17Networks_N1000.annot'));

        else
            error('yeo must be either 7 or 17.');
        end

        Yeo_ROI_Names = string(Yeo_ROI_Names)';

        %% Initialize output table

        TableYeoElectrodes = array2table(zeros(1, yeo));
        TableYeoElectrodes.Electrode_name{1} = 'label';

        TableYeoElectrodes.Properties.VariableNames = ...
            [Yeo_ROI_Names, "Electrode_name"];

        %% Load surfaces

        gL = gifti(fullfile( ...
            localDataPath, ...
            'derivatives', ...
            'freesurfer', ...
            ['sub-' sub_label], ...
            'pial.L.surf.gii'));

        gR = gifti(fullfile( ...
            localDataPath, ...
            'derivatives', ...
            'freesurfer', ...
            ['sub-' sub_label], ...
            'pial.R.surf.gii'));

        %% Process electrodes

        for ii = 1:elect_num

            LabelsYeo = strings(0, 1);

            % Select hemisphere
            if strcmpi(loc_info.hemisphere{ii}, 'L')

                g = gL;
                surface_labels.vol = verlabL;

            elseif strcmpi(loc_info.hemisphere{ii}, 'R')

                g = gR;
                surface_labels.vol = verlabR;

            end

            electrode_Yeo_areas = zeros(size(Yeo_ROI_Names, 2), 1);

            % Distance between electrode and surface vertices
            aa = sqrt(sum((g.vertices - xyz(ii,:)).^2, 2));

            % Vertices within specified distance
            near_areas = find(aa < mmDistance);

            nna = length(near_areas);

            areasNumberYeo = temp.table(2:end, 5)';

            %% Count Yeo labels around electrode

            for proba_area = 1:nna

                Yeo_area = surface_labels.vol(near_areas(proba_area));

                if Yeo_area > 0

                    for j = 1:yeo

                        if Yeo_area == areasNumberYeo(j)

                            electrode_Yeo_areas(j) = ...
                                electrode_Yeo_areas(j) + 1;

                            LabelsYeo(end + 1, 1) = Yeo_ROI_Names(j);
                        end
                    end
                end
            end

            eWa = electrode_Yeo_areas';

            njj = njj + 1;

            %% Add electrode to output table

            TableYeoElectrodes.Electrode_name{njj} = ...
                string(loc_info.name{ii});

            uniqueLabels = unique(LabelsYeo, 'stable');
            
            if isempty(uniqueLabels)
                TableYeoElectrodes.Yeo_Labels{njj} = '';
            else
                TableYeoElectrodes.Yeo_Labels{njj} = ...
                    char(strjoin(uniqueLabels, ', '));
            end

            TableYeoElectrodes{njj, 1:size(eWa,2)} = eWa;

        end

        %% Save output

        if yeo == 7
            outputDir = fullfile( ...
                localDataPath, ...
                'derivatives', ...
                'Yeo7Electrodes');

        elseif yeo == 17
            outputDir = fullfile( ...
                localDataPath, ...
                'derivatives', ...
                'Yeo17Electrodes');
        end

        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        filename = fullfile( ...
            outputDir, ...
            ['sub-' sub_label '_TableYeoElectrodes.xlsx']);

        writetable(TableYeoElectrodes, filename, 'Sheet', 1);

    end
end