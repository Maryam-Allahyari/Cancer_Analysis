"""Known marker genes (Ensembl IDs, as used in the TPM file) and their names.

The liver genes are highly expressed in normal liver. I checked that they are
much higher in liver biopsies than in other biopsies in this dataset, which is
what we would expect if some normal liver tissue ends up in the sample.
"""

LIVER_GENES = {
    "ENSG00000163631": "ALB",
    "ENSG00000118137": "APOA1",
    "ENSG00000158874": "APOA2",
    "ENSG00000091583": "APOH",
    "ENSG00000171564": "FGB",
    "ENSG00000257017": "HP",
    "ENSG00000091513": "TF",
}

# classic tissue-of-origin markers, used to label plots
TISSUE_MARKERS = {
    "ENSG00000091831": "ESR1",
    "ENSG00000107485": "GATA3",
    "ENSG00000165556": "CDX2",
    "ENSG00000136352": "NKX2-1",
    "ENSG00000168484": "SFTPB",
}

GENE_NAMES = {**LIVER_GENES, **TISSUE_MARKERS}


def pretty_name(feature):
    """Turn a feature name into something readable for plots."""
    return GENE_NAMES.get(feature, feature)
