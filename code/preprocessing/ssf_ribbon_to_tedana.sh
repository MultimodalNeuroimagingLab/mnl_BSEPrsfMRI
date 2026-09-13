#!/bin/bash

set -euo pipefail

# Resample the fMRIPrep cortical ribbon mask to the tedana BOLD space
# and create a binary gray-matter mask.
#
# Author: Maria Guadalupe Yanez Ramos
# Developed with scientific and technical guidance from Dora Hermes
# and the Multimodal Neuroimaging Lab (MNL) team.
# May 2026
#
# Usage:
#   bash ssf_ribbon_to_tedana.sh <subject-label>
#
# Example:
#   bash ssf_ribbon_to_tedana.sh 01
#
# The subject label should be provided without the "sub-" prefix.

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <subject-label>"
    exit 1
fi

SUBJ="${1#sub-}"

# Determine project location from this script:
# <project>/code/preprocessing/ssf_ribbon_to_tedana.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

DERIVATIVES_DIR="${PROJECT_DIR}/data/derivatives"

ANAT_DIR="${DERIVATIVES_DIR}/fmriprep/sub-${SUBJ}/ses-compact3T01/anat"
TEDANA_DIR="${DERIVATIVES_DIR}/tedana/sub-${SUBJ}"

IN_MASK="${ANAT_DIR}/sub-${SUBJ}_ses-compact3T01_desc-ribbon_mask.nii.gz"

REF_BOLD="${TEDANA_DIR}/sub-${SUBJ}_ses-compact3T01_task-rest_run-01_desc-optcomAccepted_bold.nii.gz"

OUT_MASK="${ANAT_DIR}/sub-${SUBJ}_ses-compact3T01_desc-ribbon_mask_tedanaSpace.nii.gz"

OUT_BIN="${ANAT_DIR}/sub-${SUBJ}_ses-compact3T01_desc-ribbon_mask_tedanaSpace_bin.nii.gz"


# Check required commands.
command -v flirt >/dev/null 2>&1 || {
    echo "Error: FSL flirt was not found on the PATH."
    exit 1
}

command -v fslmaths >/dev/null 2>&1 || {
    echo "Error: FSL fslmaths was not found on the PATH."
    exit 1
}


# Check required input files.
if [ ! -f "${IN_MASK}" ]; then
    echo "Error: ribbon mask not found:"
    echo "${IN_MASK}"
    exit 1
fi

if [ ! -f "${REF_BOLD}" ]; then
    echo "Error: tedana BOLD file not found:"
    echo "${REF_BOLD}"
    exit 1
fi


echo "Processing sub-${SUBJ}"
echo "Resampling cortical ribbon mask to tedana space..."

flirt \
    -in "${IN_MASK}" \
    -ref "${REF_BOLD}" \
    -out "${OUT_MASK}" \
    -applyxfm \
    -usesqform


echo "Creating binary ribbon mask..."

fslmaths \
    "${OUT_MASK}" \
    -thr 0.5 \
    -bin \
    "${OUT_BIN}"


echo "Done."
echo "Resampled mask:"
echo "${OUT_MASK}"
echo "Binary mask:"
echo "${OUT_BIN}"
