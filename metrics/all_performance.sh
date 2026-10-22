#!/bin/bash

# Change it manually:
IN_FILE="adata_jakob_raw_log1pNorm_hvg_3000"   # without .h5ad
META_DIR="/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/BASS/predicted_clusters/3000_hvg_log1pNorm_20PC_R3_C_celltype_r5_predclusters/" # change directory
OUT_DIR="/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/BASS/" # change method name

# Automatic:
FOLDER_NAME="performance_$(basename "${META_DIR%/}")"
JOB_DIR="${OUT_DIR}/${FOLDER_NAME}"

echo "Creating new directory: $JOB_DIR"
mkdir -p "$JOB_DIR"

shopt -s nullglob
csv_files=("$META_DIR"/*.csv)

if [ ${#csv_files[@]} -eq 0 ]; then
    echo "Error: .csv files not found in $META_DIR"
    exit 1
fi

# Submiting sbatch jobs: CHANGE INPUTS MANUALLY

for file_path in "${csv_files[@]}"; do
    file_name=$(basename "$file_path" .csv)

    sbatch --job-name="performance_${file_name}" \
           --account=limitedAccessBladder \
           --time=00:30:00 \
           --mem=30GB \
           --ntasks=1 \
           --cpus-per-task=8 \
           --wrap="python get_performance_metrics.py \
                   --in_file $IN_FILE \
                   --out_path $JOB_DIR \
                   --meta_dir $META_DIR \
                   --file_predicted_clustering $file_name \
                   --method BASS \
                   --spatial_key spatial_banksy \
                   --annotation_key celltype_r5 \
                   --spatial_x CenterX_global_px \
                   --spatial_y CenterY_global_px"
done

echo "Jobs have been sent. All in: $JOB_DIR"