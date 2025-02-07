#!/bin/bash

OUT_DIR=/path/to/results
SOURCE_SCORE_DIR=/path/to/PRS/data
SAMPLE_SHEET=/path/to/PGS_samplesheet.csv
ANCESTRY_FILE=/path/to/1000genomes/ancestry/file/pgsc_1000G_v1.tar.zst
target_build="38"

mkdir -p $OUT_DIR

echo $SAMPLE_SHEET

nextflow run pgscatalog/pgsc_calc \
    -profile singularity \
    --input $SAMPLE_SHEET \
    --scorefile "${SOURCE_SCORE_DIR}/*_hg${target_build}.txt" \
    --target_build "GRCh${target_build}" \
    --outdir ${OUT_DIR} \
    --max_cpus 16 \
    --max_memory 16.GB \
    --run_ancestry $ANCESTRY_FILE \
    --min_overlap 0 | tee ${OUT_DIR}/pipeline.log 2>&1 
   

nextflow run pgscatalog/pgsc_calc \
    -profile singularity \
    --input $SAMPLE_SHEET \
    --pgs_id PGS003725,PGS000507\
    --target_build "GRCh${target_build}" \
    --outdir ${OUT_DIR} \
    --max_cpus 16 \
    --max_memory 16.GB \
    --run_ancestry $ANCESTRY_FILE \
    --min_overlap 0 | tee ${OUT_DIR}/pipeline.log 2>&1 