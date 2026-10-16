#!/usr/bin/env bash
# Run every step one after the other without Nextflow (handy for debugging).
# Assumes the Python packages from requirements.txt are installed.
set -euo pipefail

RAW=Data/Raw
PROC=Data/Processed
RES=results

python src/clean.py \
    --expression $RAW/POG570_TPM_expression.txt.gz \
    --mutations $RAW/POG570_small_mutations.txt.gz \
    --clinical $RAW/Table_S1_Demographics.xlsx \
    --out-dir $PROC

python src/features.py \
    --mutations $PROC/mutations_coding.csv --tmb $PROC/tmb.csv \
    --clinical $PROC/clinical.csv --genes resources/driver_genes.txt \
    --out $PROC/mutation_features.csv

python src/explore.py \
    --expression $PROC/expression.parquet --mutation-features $PROC/mutation_features.csv \
    --clinical $PROC/clinical.csv --out-dir $RES/explore

python src/train.py \
    --expression $PROC/expression.parquet --mutation-features $PROC/mutation_features.csv \
    --clinical $PROC/clinical.csv --out-dir $RES/models

PYTHONPATH=src python src/explain.py \
    --model $RES/models/final_model.joblib \
    --expression $PROC/expression.parquet --mutation-features $PROC/mutation_features.csv \
    --clinical $PROC/clinical.csv \
    --test-predictions $RES/models/test_predictions.csv \
    --oof-predictions $RES/models/oof_predictions.csv \
    --out-dir $RES/explain
