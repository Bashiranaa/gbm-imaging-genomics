#!/usr/bin/env bash

RAW_ROOT=~/UPENN-GBM
NII_ROOT=~/nifti_converted
mkdir -p "$NII_ROOT"

for subj_dir in "${RAW_ROOT}"/UPENN-GBM-*; do
  subjID=$(basename "${subj_dir}")
  out_dir="${NII_ROOT}/${subjID}"
  mkdir -p "${out_dir}"

  echo "Converting ${subjID}..."

  dcm2niix -z y \
           -f "${subjID}_T1pre" \
           -o "${out_dir}" \
           "${subj_dir}" \
  && echo "Done with ${subjID}" \
  || echo "Problem converting ${subjID}"
done

