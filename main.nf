#!/usr/bin/env nextflow

/*
 * Cancer origin prediction pipeline
 * CLEAN -> FEATURES -> EXPLORE / TRAIN -> EXPLAIN
 *
 * Run:  nextflow run main.nf -profile docker
 */

nextflow.enable.dsl = 2

// the Python code is passed to every process, so it also works inside containers
src_dir = file("${projectDir}/src")


process CLEAN {
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path src
    path expression
    path mutations
    path clinical

    output:
    path "clinical.csv",          emit: clinical
    path "expression.parquet",    emit: expression
    path "mutations_coding.csv",  emit: coding
    path "tmb.csv",               emit: tmb

    script:
    """
    python ${src}/clean.py \\
        --expression ${expression} \\
        --mutations ${mutations} \\
        --clinical ${clinical} \\
        --out-dir . \\
        --min-patients ${params.min_patients} \\
        --min-af ${params.min_af} \\
        --min-alt-reads ${params.min_alt_reads}
    """
}


process FEATURES {
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path src
    path coding
    path tmb
    path clinical
    path genes

    output:
    path "mutation_features.csv"

    script:
    """
    python ${src}/features.py \\
        --mutations ${coding} --tmb ${tmb} --clinical ${clinical} \\
        --genes ${genes} --out mutation_features.csv
    """
}


process EXPLORE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path src
    path expression
    path mut_features
    path clinical

    output:
    path "explore"

    script:
    """
    python ${src}/explore.py \\
        --expression ${expression} --mutation-features ${mut_features} \\
        --clinical ${clinical} --out-dir explore
    """
}


process TRAIN {
    publishDir "${params.outdir}", mode: 'copy'
    cpus 4

    input:
    path src
    path expression
    path mut_features
    path clinical

    output:
    path "models",                               emit: all
    path "models/final_model.joblib",            emit: model
    path "models/test_predictions.csv",          emit: test_pred
    path "models/oof_predictions.csv",           emit: oof_pred

    script:
    """
    python ${src}/train.py \\
        --expression ${expression} --mutation-features ${mut_features} \\
        --clinical ${clinical} --out-dir models \\
        --n-genes ${params.n_genes} --seed ${params.seed} --tune-iter ${params.tune_iter}
    """
}


process EXPLAIN {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path src
    path model
    path expression
    path mut_features
    path clinical
    path test_pred
    path oof_pred

    output:
    path "explain"

    script:
    """
    export PYTHONPATH=${src}
    python ${src}/explain.py \\
        --model ${model} --expression ${expression} \\
        --mutation-features ${mut_features} --clinical ${clinical} \\
        --test-predictions ${test_pred} --oof-predictions ${oof_pred} \\
        --out-dir explain
    """
}


workflow {
    expression = file("${params.raw_dir}/POG570_TPM_expression.txt.gz", checkIfExists: true)
    mutations  = file("${params.raw_dir}/POG570_small_mutations.txt.gz", checkIfExists: true)
    clinical   = file("${params.raw_dir}/Table_S1_Demographics.xlsx", checkIfExists: true)
    genes      = file(params.driver_genes, checkIfExists: true)

    CLEAN(src_dir, expression, mutations, clinical)
    FEATURES(src_dir, CLEAN.out.coding, CLEAN.out.tmb, CLEAN.out.clinical, genes)

    EXPLORE(src_dir, CLEAN.out.expression, FEATURES.out, CLEAN.out.clinical)
    TRAIN(src_dir, CLEAN.out.expression, FEATURES.out, CLEAN.out.clinical)

    EXPLAIN(src_dir, TRAIN.out.model, CLEAN.out.expression, FEATURES.out,
            CLEAN.out.clinical, TRAIN.out.test_pred, TRAIN.out.oof_pred)
}
