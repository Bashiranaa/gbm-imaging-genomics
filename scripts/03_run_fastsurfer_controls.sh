#!/bin/bash

INPUT_DIR="/home/student24/nifti_reoriented_controls"
OUTPUT_DIR="/home/student24/Desikan_Segmentation_controls"

mkdir -p "$OUTPUT_DIR"

# Loop over each subject folder named sub-XXX inside INPUT_DIR
for sample in "$INPUT_DIR"/sub-*; do
    # basename "$sample" will be something like "sub-001"
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

