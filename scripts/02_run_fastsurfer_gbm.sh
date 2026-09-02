#!/bin/bash

INPUT_DIR="/home/student24/nifti_reoriented"
OUTPUT_DIR="/home/student24/Desikan_Segmentation"

mkdir -p "$OUTPUT_DIR"

for sample in "$INPUT_DIR"/UPENN-GBM-*; do
    sample_id=$(basename "$sample")
    echo "Processing $sample_id..."

    docker run --rm \
  -u $(id -u):$(id -g) \
  -v "$sample":/input \
  -v "$OUTPUT_DIR":/output \
  -v /home/student24/Downloads/license.txt:/license.txt \
  deepmi/fastsurfer \
  --t1 /input/first_volume.nii.gz \
  --sid "$sample_id" \
  --sd /output \
  --fs_license /license.txt \
  --threads 2
done

