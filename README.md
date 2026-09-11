# mnl_BSEPrsfMRI
integrating resting state fMRI with BSEPS

# BSEP–rs-fMRI Analysis Workflow

This repository contains the analysis workflow used to integrate brain stimulation evoked potentials (BSEPs) with resting-state fMRI functional connectivity.

The workflow includes:

- rs-fMRI preprocessing with fMRIPrep and tedana
- subject-specific Yeo 7-network projection
- electrode registration to rs-fMRI anatomical space
- stimulation-pair selection
- BSEP/CRP preprocessing
- rs-fMRI ROI time-series extraction
- integration of electrophysiology and rs-fMRI data
- functional connectivity analysis
- within-subject FC–EC correlation visualization

---

# PART 1 — Repository setup and required data

## 1. Clone the repository

Open a terminal and choose the directory where the project will be stored.

```bash
cd /path/to/projects

git clone https://github.com/MultimodalNeuroimagingLab/mnl_BSEPrsfMRI.git

cd mnl_BSEPrsfMRI
```

The repository should contain approximately:

```text
mnl_BSEPrsfMRI/
├── code/
│   ├── functions/
│   ├── preprocessing/
│   ├── ssf_01_electrode_network_mapping.m
│   ├── ssf_02_subject_selection.m
│   ├── ssf_03_preprocess_bsep_crp.m
│   ├── ssf_04_extract_rsfmri_roi_timeseries.m
│   ├── ssf_05_integrate_bsep_rsfmri_data.m
│   ├── ssf_06_compute_fc_analysis_table.m
│   └── ssf_07_plot_within_subject_fc_ec_correlations.m
├── LICENSE
└── README.md
```

The `data/` directory is not stored in GitHub.

---

## 2. Create the local data directory

From the repository root:

```bash
mkdir -p data
```

All subject-specific data and generated derivatives should remain inside:

```text
mnl_BSEPrsfMRI/data/
```

---

## 3. Copy the required subject data

For each subject, copy the BIDS iEEG directory required for the BSEP/CCEP analysis.

Expected structure:

```text
data/
└── sub-<label>/
    └── ses-ieeg01/
        └── ieeg/
            ├── sub-<label>_ses-ieeg01_electrodes.tsv
            ├── sub-<label>_ses-ieeg01_task-ccep_run-*_events.tsv
            ├── sub-<label>_ses-ieeg01_task-ccep_run-*_channels.tsv
            └── sub-<label>_ses-ieeg01_task-ccep_run-*_ieeg.mefd/
```

Copy all candidate CCEP runs rather than only the runs expected to be selected.

Stage 02 examines the available stimulation data and determines which stimulation pairs and runs are used downstream.

---

## 4. Required MRI and electrode-localization derivatives

Before the MATLAB workflow is started, the following subject-specific inputs should be available.

Recommended structure:

```text
data/
└── derivatives/
    ├── T1electrodeCoordinates/
    │   └── sub-<label>/
    │       └── sub-<label>_ses-compact3T01_desc-acpc_T1w.nii
    │
    ├── fmriprep/
    │   └── sub-<label>/
    │       └── ses-compact3T01/
    │           ├── anat/
    │           │   ├── sub-<label>_ses-compact3T01_desc-preproc_T1w.nii.gz
    │           │   └── sub-<label>_ses-compact3T01_desc-ribbon_mask_tedanaSpace_bin.nii.gz
    │           └── func/
    │               └── multi-echo preprocessed BOLD files
    │
    ├── freesurfer/
    │   └── sub-<label>/
    │       └── ...
    │
    └── tedana/
        └── sub-<label>/
            └── sub-<label>_ses-compact3T01_task-rest_run-01_desc-optcomAccepted_bold.nii.gz
```

`T1electrodeCoordinates` contains the anatomical T1 image in the same space as the original electrode coordinates.

The FreeSurfer subject directory must contain the subject-specific Yeo 7-network annotation files:

```text
label/lh.Yeo2011_7Networks_N1000.annot
label/rh.Yeo2011_7Networks_N1000.annot
```

---

# PART 2 — rs-fMRI preprocessing required before MATLAB

The following steps must be completed for each subject before running MATLAB Stages 01–07.

---

## 1. Run fMRIPrep

The resting-state fMRI data are processed with fMRIPrep.

For multi-echo data, it is important to include:

```text
--me-output-echos
```

because the individual preprocessed echoes are later used as inputs to tedana.

An example fMRIPrep command is provided in:

```text
code/preprocessing/fmrprepServer.txt
```


---

## 2. Remove the first 5 volumes from each preprocessed echo

Before running tedana, remove the first five volumes from each preprocessed echo.

Example for one echo:

```bash
fslroi \
    sub-01_ses-compact3T01_task-rest_run-01_echo-01_desc-preproc_bold.nii.gz \
    temp_echo-01.nii.gz \
    5 -1

mv \
    temp_echo-01.nii.gz \
    sub-01_ses-compact3T01_task-rest_run-01_echo-01_desc-preproc_bold.nii.gz
```

Repeat this process for all echoes used by tedana.

For example:

```text
echo-01
echo-02
echo-03
```

Important: `Auto_tedana_last.py` searches for files ending in:

```text
_desc-preproc_bold.nii.gz
```

Therefore, when tedana is run, the five-volume-trimmed images must retain this filename pattern.

---

## 3. Run tedana

The repository contains:

```text
code/preprocessing/Auto_tedana_last.py
```

The script:
- identifies the preprocessed multi-echo BOLD files
- identifies the corresponding BIDS JSON files
- reads the echo times from the JSON files
- identifies the fMRIPrep brain mask
- creates the subject tedana output directory
- runs `tedana_workflow`

```bash
PROJECT=/path/to/mnl_BSEPrsfMRI

python3 \
    "$PROJECT/code/preprocessing/Auto_tedana_last.py" \
    --fmriprepDir "$PROJECT/data/derivatives/fmriprep" \
    --bidsDir "$PROJECT/data" \
    --cores 4 \
    --subjectName 01
```

Verify that the final tedana BOLD image required by the MATLAB workflow exists:

```text
data/derivatives/tedana/sub-01/
sub-01_ses-compact3T01_task-rest_run-01_desc-optcomAccepted_bold.nii.gz
```

---

## 4. Transform the cortical ribbon mask to tedana space

After tedana, the cortical gray-matter ribbon mask must be resampled to the tedana BOLD space and binarized.

The repository includes:

```text
code/preprocessing/ribbon2tedana.sh
```

The main operations are:

```bash
flirt \
    -in <ribbon_mask> \
    -ref <tedana_bold> \
    -out <ribbon_tedana_space> \
    -applyxfm \
    -usesqform

fslmaths \
    <ribbon_tedana_space> \
    -thr 0.5 \
    -bin \
    <ribbon_tedana_space_bin>
```

The final file required by MATLAB Stage 04 is:

```text
data/derivatives/fmriprep/sub-<label>/ses-compact3T01/anat/
sub-<label>_ses-compact3T01_desc-ribbon_mask_tedanaSpace_bin.nii.gz
```

This mask is used to restrict electrode ROIs to gray-matter voxels.

### Important

`ribbon2tedana.sh` may contain local path definitions from the original development environment.

---

## 5. Project the Yeo 7-network atlas to the subject surface

Use FreeSurfer `mri_surf2surf` to project the Yeo 7-network atlas from `fsaverage` to the individual subject surface.

The repository includes:

```text
code/preprocessing/Yeo7_surf2surf.sh
```

The script performs the projection separately for the left and right hemispheres.

Example:

```bash
bash code/preprocessing/Yeo7_surf2surf.sh 01
```

Verify that the subject FreeSurfer directory contains:

```text
label/lh.Yeo2011_7Networks_N1000.annot
label/rh.Yeo2011_7Networks_N1000.annot
```

These subject-specific annotations are required before the electrode network assignment performed in MATLAB Stage 01.

### Important


Before running `Yeo7_surf2surf.sh` on a new computer, inspect and update:

```text
SUBJECTS_DIR
fsaverage location
subject FreeSurfer directory
```

as needed.

---

## 6. Verify required inputs before starting MATLAB

Before running Stage 01, verify that each subject has the required inputs.

At minimum:

```text
data/
├── sub-<label>/
│   └── ses-ieeg01/
│       └── ieeg/
│           ├── electrodes.tsv
│           ├── events.tsv
│           ├── channels.tsv
│           └── MEF data
│
└── derivatives/
    ├── T1electrodeCoordinates/
    │   └── sub-<label>/
    │       └── sub-<label>_ses-compact3T01_desc-acpc_T1w.nii
    │
    ├── fmriprep/
    │   └── sub-<label>/
    │       └── ses-compact3T01/
    │           └── anat/
    │               ├── sub-<label>_ses-compact3T01_desc-preproc_T1w.nii.gz
    │               └── sub-<label>_ses-compact3T01_desc-ribbon_mask_tedanaSpace_bin.nii.gz
    │
    ├── freesurfer/
    │   └── sub-<label>/
    │       └── label/
    │           ├── lh.Yeo2011_7Networks_N1000.annot
    │           └── rh.Yeo2011_7Networks_N1000.annot
    │
    └── tedana/
        └── sub-<label>/
            └── sub-<label>_ses-compact3T01_task-rest_run-01_desc-optcomAccepted_bold.nii.gz
```

Once these files are available, the MATLAB workflow can be started.

---

# PART 3 — MATLAB setup and BSEP–rs-fMRI workflow

## 1. Install MATLAB dependencies

The workflow requires MATLAB together with external neuroimaging and iEEG dependencies.

Required or currently used dependencies include:

```text
Vistasoft
SPM
mnl_ieegBasics
matmef
CRP analysis code
MATLAB Signal Processing Toolbox
MATLAB Statistics and Machine Learning Toolbox
```

Additional FreeSurfer and FSL command-line tools are required for the preprocessing steps described above.

---

## 2. Add the required code to the MATLAB path

Start MATLAB and go to the cloned repository:

```matlab
cd('/path/to/mnl_BSEPrsfMRI')
```

Add the repository code:

```matlab
projectDir = pwd;

addpath(genpath(fullfile(projectDir, 'code')));
```

Then add the locally installed external dependencies.

Example:

```matlab
addpath(genpath('/path/to/vistasoft'));
addpath(genpath('/path/to/mnl_ieegBasics'));
addpath(genpath('/path/to/matmef'));
addpath('/path/to/spm12');
addpath(genpath('/path/to/crp_scripts'));
```

The exact paths depend on the local installation.

---

## 3. Verify MATLAB dependencies

Before running the workflow, verify that MATLAB can locate the required functions.

```matlab
which niftiRead
which niftiWrite

which spm_coreg
which spm_matrix

which read_annotation

which ieeg_FSsurf2T1space
which ieeg_eldist2pial2white

which ccep_PreprocessMef
which ssf_CAR64blocks_percent
which CRP_method

which fitlm
```

Each command should return a valid file path.

---

# PART 4 — Run MATLAB Stages 01–07

Run the stages sequentially.

---

## Stage 01 — Electrode network mapping

Run:

```matlab
ssf_01_electrode_network_mapping
```

Stage 01 automatically discovers subjects from:

```text
data/derivatives/T1electrodeCoordinates/sub-*
```

For each subject, the stage:

- registers electrode coordinates to rs-fMRI T1 space
- creates T1-space cortical surface representations
- assigns Yeo 7-network labels to electrodes
- generates cortical-network visualizations
- computes relative gray/white-matter distance

Main outputs:

```text
data/derivatives/electrodes2rsfMRI/
data/derivatives/Yeo7Electrodes/
data/derivatives/gm_wmDistance/
```

Example transformed electrode file:

```text
data/derivatives/electrodes2rsfMRI/sub-<label>/
sub-<label>_ses-ieeg01_space-rsT1_electrodes.tsv
```

Inspect the registration and electrode/network localization quality-control outputs before continuing.

---

## Stage 02 — Subject and stimulation-pair selection

Run:

```matlab
ssf_02_subject_selection
```

Stage 02 automatically discovers subjects with Stage 01 electrode derivatives and identifies eligible bipolar stimulation pairs.

Selection criteria include:

- available CCEP/BSEP stimulation data
- Yeo Visual or Control network membership
- gray/white-matter relative distance
- seizure-onset-zone exclusion
- availability of processable CCEP runs

Main output:

```text
data/vcm_Subjects.mat
```

Additional subject-level information is stored in:

```text
data/derivatives/loc_info/
└── sub-<label>/
    └── loc_info.mat
```

---

## Stage 03 — BSEP and CRP preprocessing

Run:

```matlab
ssf_03_preprocess_bsep_crp
```

For each stimulation pair selected in Stage 02, this stage reads the corresponding:

```text
_events.tsv
_channels.tsv
_ieeg.mefd/
```

The stage performs:

- BSEP/CCEP trial loading
- baseline correction
- recording-channel selection
- common-average rereferencing
- notch filtering
- CRP analysis

Recording-channel selection uses information from `channels.tsv` and electrode localization information, including SOZ labels and distance from stimulation contacts.

Main outputs:

```text
data/derivatives/whimstim/
data/derivatives/preprocPlots/
data/derivatives/crpMPlots/
```

Review the preprocessing and CRP QC figures before continuing.

---

## Stage 04 — rs-fMRI electrode ROI time series

Run:

```matlab
ssf_04_extract_rsfmri_roi_timeseries
```

This stage uses:

```text
data/vcm_Subjects.mat
data/derivatives/electrodes2rsfMRI/
data/derivatives/tedana/
data/derivatives/fmriprep/
```

It creates small gray-matter-constrained ROIs around individual electrode contacts and stimulation pairs and extracts resting-state fMRI time series.

Main outputs:

```text
data/derivatives/maskROIS/
└── sub-<label>/
    └── ROI_*.nii
```

and:

```text
data/derivatives/timeSeries/
└── sub-<label>/
    └── timeSeries.mat
```

`timeSeries.mat` contains the individual-contact and stimulation-pair ROI time-series structures used in Stage 05.

---

## Stage 05 — Integrate BSEP and rs-fMRI data

Run:

```matlab
ssf_05_integrate_bsep_rsfmri_data
```

This stage combines:

- BSEP/CRP results from Stage 03
- rs-fMRI time series from Stage 04
- electrode localization information
- Yeo network information
- gray/white-matter distance
- stimulation-pair information

Main output:

```text
data/derivatives/vcm_data/
└── sub-<label>/
    └── <label>_<stimulation-pair>.mat
```

These files contain the integrated per-recording-channel information used for the functional-connectivity analysis.

---

## Stage 06 — Functional-connectivity analysis table

Run:

```matlab
ssf_06_compute_fc_analysis_table
```

Stage 06 calculates resting-state functional connectivity between each stimulation seed and recording-contact ROI.

It combines functional connectivity with BSEP/CRP measurements and generates a long-format analysis table.

Main output:

```text
data/allRows.mat
```

The output table includes information such as:

```text
subject
stim_pair
stim_network
rec_channel
rec_network
insideStimNetwork
gm_wm_relativeDistance
Add2Analyses
rPearson
pPearson
zFisher
median_cod
mean_cod_noOutliers
sigRes
```

Stage 06 also assigns independent Yeo-network membership indicators for recording contacts and determines whether each recording contact lies inside or outside the stimulated network.

Inspect the summary printed in MATLAB after completion.

---

## Stage 07 — Within-subject FC–EC correlations

Run:

```matlab
ssf_07_plot_within_subject_fc_ec_correlations
```

Stage 07 examines the relationship between:

```text
FC = resting-state functional connectivity, Fisher z
EC = BSEP median coefficient of determination (CoD)
```

The visualization distinguishes same-network recording contacts from the complete set of valid recording contacts.

Main figure output:

```text
data/derivatives/figures/fc_ec_correlations/
└── ssf_within_subject_fc_ec_correlations.png
```

---

# PART 5 — Outputs from a clean rerun

For a clean reproducibility run, downstream analysis derivatives should be generated by Stages 01–07 rather than copied from a previous analysis.

Do not copy previously generated versions of:

```text
data/vcm_Subjects.mat
data/allRows.mat

data/derivatives/electrodes2rsfMRI/
data/derivatives/Yeo7Electrodes/
data/derivatives/gm_wmDistance/
data/derivatives/loc_info/
data/derivatives/whimstim/
data/derivatives/preprocPlots/
data/derivatives/crpMPlots/
data/derivatives/maskROIS/
data/derivatives/timeSeries/
data/derivatives/vcm_data/
data/derivatives/figures/
```

After completing Stage 07, the important generated products should include:

```text
data/
├── vcm_Subjects.mat
├── allRows.mat
└── derivatives/
    ├── electrodes2rsfMRI/
    ├── Yeo7Electrodes/
    ├── gm_wmDistance/
    ├── loc_info/
    ├── whimstim/
    ├── preprocPlots/
    ├── crpMPlots/
    ├── maskROIS/
    ├── timeSeries/
    ├── vcm_data/
    └── figures/
```

## MATLAB and computing environment

The analysis pipeline was tested using MATLAB R2024a on Apple silicon macOS systems, including:

MacBook Pro, Apple M4 Max, 128 GB RAM
Mac Studio, Apple M3 Ultra, 96 GB RAM

Testing was performed under macOS Tahoe 26.6.2.

Because several cortical rendering and connectivity-figure steps can require substantial memory, systems with sufficient RAM are recommended.
