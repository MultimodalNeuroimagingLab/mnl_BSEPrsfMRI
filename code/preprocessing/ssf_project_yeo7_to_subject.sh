#!/bin/bash

set -euo pipefail

# Project the Yeo 7-network cortical atlas from FreeSurfer fsaverage
# to an individual subject surface.
#
# Author: Maria Guadalupe Yanez Ramos
# Developed with scientific and technical guidance from Dora Hermes
# and the Multimodal Neuroimaging Lab (MNL) team.
# May 2026
#
# Usage:
#   bash ssf_project_yeo7_to_subject.sh <subject-label>
#
# Example:
#   bash ssf_project_yeo7_to_subject.sh 01
#
# The subject label should be provided without the "sub-" prefix.

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <subject-label>"
    exit 1
fi

SUBJ="${1#sub-}"

# Determine project location from this script:
# <project>/code/preprocessing/ssf_project_yeo7_to_subject.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

export SUBJECTS_DIR="${PROJECT_DIR}/data/derivatives/freesurfer"

SOURCE_SUBJECT="fsaverage"
TARGET_SUBJECT="sub-${SUBJ}"

SOURCE_LABEL_DIR="${SUBJECTS_DIR}/${SOURCE_SUBJECT}/label"
TARGET_LABEL_DIR="${SUBJECTS_DIR}/${TARGET_SUBJECT}/label"


# Check FreeSurfer command.
command -v mri_surf2surf >/dev/null 2>&1 || {
    echo "Error: FreeSurfer mri_surf2surf was not found on the PATH."
    exit 1
}


# Check required directories.
if [ ! -d "${SUBJECTS_DIR}/${SOURCE_SUBJECT}" ]; then
    echo "Error: fsaverage directory not found:"
    echo "${SUBJECTS_DIR}/${SOURCE_SUBJECT}"
    exit 1
fi

if [ ! -d "${SUBJECTS_DIR}/${TARGET_SUBJECT}" ]; then
    echo "Error: subject FreeSurfer directory not found:"
    echo "${SUBJECTS_DIR}/${TARGET_SUBJECT}"
    exit 1
fi

mkdir -p "${TARGET_LABEL_DIR}"


# Project Yeo atlas for both hemispheres.
for HEMI in lh rh; do

    SOURCE_ANNOT="${SOURCE_LABEL_DIR}/${HEMI}.Yeo2011_7Networks_N1000.annot"

    TARGET_ANNOT="${TARGET_LABEL_DIR}/${HEMI}.Yeo2011_7Networks_N1000.annot"

    if [ ! -f "${SOURCE_ANNOT}" ]; then
        echo "Error: Yeo annotation not found:"
        echo "${SOURCE_ANNOT}"
        exit 1
    fi

    echo "Projecting Yeo 7-network atlas for ${TARGET_SUBJECT}, ${HEMI}..."

    mri_surf2surf \
        --srcsubject "${SOURCE_SUBJECT}" \
        --trgsubject "${TARGET_SUBJECT}" \
        --hemi "${HEMI}" \
        --sval-annot "${SOURCE_ANNOT}" \
        --tval "${TARGET_ANNOT}"

done


echo "Done."
echo "Generated annotations:"
echo "${TARGET_LABEL_DIR}/lh.Yeo2011_7Networks_N1000.annot"
echo "${TARGET_LABEL_DIR}/rh.Yeo2011_7Networks_N1000.annot"
