# Predicting Cancer Origin from Tumour Genomics

Can a machine learning model tell where a metastatic tumour started, using only the tumour's gene expression and mutations?

This project tries to answer that with real patient data from the Personalized OncoGenomics (POG) program at BC Cancer (570 patients with advanced cancer, Pleasance et al., 2020). It is a work in progress. First results are below.

## Why this question matters

Treatment for cancer is usually chosen by the organ it came from. Breast cancer that has spread to the liver is still treated as breast cancer. In about 3–5% of patients, though, doctors find metastases but cannot identify the primary site, even after a full work-up. This is called cancer of unknown primary (CUP), and it usually has a poor outlook (Pavlidis & Pentheroudakis, 2012).

Tumours tend to keep some molecular "memory" of their tissue of origin. Several groups have used this to build classifiers from RNA-seq data (Grewal et al., 2019), from targeted DNA panels (Penson et al., 2020; Moon et al., 2023) and from whole-genome mutation patterns (Jiao et al., 2020). This project is a smaller, fully reproducible version of the same idea that I can run end to end on a laptop.

In the POG cohort the true cancer type of every patient is known, so this is not a CUP study. The point is to see how well the origin can be recovered from metastatic biopsies when we already know the right answer.

## Questions

1. How accurately can gene expression predict the cancer type of a metastatic biopsy?
2. Do mutation features (driver gene mutations, tumour mutation burden) add anything on top of expression, or do they work on their own?
3. Does the biopsy site get in the way? A sample taken from a liver metastasis contains some normal liver tissue, and a model could learn "liver" instead of the real primary.
4. Which genes drive the predictions, and do they match known tissue markers?

## Data

All data come from the POG570 release by BC Cancer's Genome Sciences Centre:
<https://www.bcgsc.ca/downloads/POG570/>

| File | Used for |
|---|---|
| `POG570_TPM_expression.txt.gz` | Gene expression features |
| `POG570_small_mutations.txt.gz` | Mutation features (somatic SNVs and indels, hg19) |
| `Table_S1_Demographics.xlsx` | Cancer type (the label) and patient details |

**The data are not included in this repository.** The release is for research use only and may not be redistributed. Download the files yourself and put them in `Data/Raw/`. See [docs/data.md](docs/data.md) for file details, cleaning rules and the terms of use.

## Approach

The plan, step by step:

1. **Load and clean.** Check for duplicate patients, missing labels and impossible values. Log how many patients and genes are kept at each step.
2. **Expression features.** Transform TPM values with log2(TPM + 1) and keep the most variable genes.
3. **Mutation features.** Keep protein-changing mutations, build a patient × gene table (mutated or not) for well-known cancer driver genes, for example from the COSMIC Cancer Gene Census (Sondka et al., 2018), and add tumour mutation burden.
4. **Exploration.** PCA of the expression data coloured by cancer type and by biopsy site, plus basic statistics (e.g. does mutation burden differ between cancer types?).
5. **Models.** A multinomial logistic regression baseline, then XGBoost (Chen & Guestrin, 2016), both built with scikit-learn pipelines (Pedregosa et al., 2011). Each model is trained three ways: expression only, mutations only, and both together.
6. **Evaluation.** Stratified 5-fold cross-validation on the training set, plus a held-out test set that is only used once at the end. I report balanced accuracy, because some cancer types have many more patients than others, along with a confusion matrix. Scaling and feature selection are fitted inside the cross-validation folds so no information leaks from the test data.
7. **Explanation.** SHAP values (Lundberg & Lee, 2017) to see which genes push a sample towards each cancer type, and whether biopsy site shows up in the explanations.
8. **Reproducibility.** The whole analysis will run as a Nextflow pipeline (Di Tommaso et al., 2017) inside a Docker container, with unit tests run by GitHub Actions on every push.

## Repository layout

```
Cancer_Analysis/
├── README.md
├── main.nf                # Nextflow pipeline
├── nextflow.config        # parameters and profiles (standard, docker, slurm)
├── Dockerfile             # container with all Python packages
├── requirements.txt       # pinned package versions
├── environment.yml        # conda alternative
├── docs/
│   └── data.md            # data description, cleaning rules, terms of use
├── Data/
│   └── Raw/               # downloaded POG570 files (not committed)
├── resources/
│   └── driver_genes.txt   # driver genes used as mutation features
├── scripts/
│   ├── download_data.sh   # fetch the POG570 files
│   └── run_local.sh       # run all steps without Nextflow
├── src/
│   ├── load.py            # read expression, mutations and demographics
│   ├── clean.py           # quality checks and filtering
│   ├── features.py        # mutation features (driver genes, TMB)
│   ├── explore.py         # PCA, liver signal, statistical tests
│   ├── train.py           # models, cross-validation, test set
│   ├── explain.py         # SHAP and biopsy-site check
│   └── markers.py         # known marker genes
├── tests/                 # pytest unit tests on small fake datasets
└── .github/workflows/     # runs the tests on every push
```

## How to run

1. Download the data (about 300 MB):
   ```bash
   bash scripts/download_data.sh
   ```
2. Run the pipeline with Nextflow and Docker (needs Linux or WSL, Java 17+, Nextflow and Docker):
   ```bash
   docker build -t cancer-origin:1.0 .
   nextflow run main.nf -profile docker
   ```
   Or, without Nextflow, install the packages and run every step in order:
   ```bash
   pip install -r requirements.txt
   bash scripts/run_local.sh
   ```
3. Run the tests (these use small fake data, no download needed):
   ```bash
   pytest
   ```

Results are written to `results/`.

## First results

From the first full run (434 patients, 7 cancer types; seed 42). These numbers will change a bit if the split or settings change.

| Model | Features | CV balanced accuracy (5-fold, mean ± SD) | Test balanced accuracy (95% bootstrap CI) |
|---|---|---|---|
| Logistic regression | Expression | 0.88 ± 0.05 | 0.89 (0.79–0.97) |
| Logistic regression | Mutations | 0.46 ± 0.05 | 0.51 (0.42–0.61) |
| Logistic regression | Both | 0.87 ± 0.05 | 0.89 (0.79–0.97) |
| XGBoost | Expression | 0.83 ± 0.05 | 0.85 (0.75–0.95) |
| XGBoost | Mutations | 0.43 ± 0.02 | 0.46 (0.37–0.57) |
| XGBoost | Both | 0.83 ± 0.03 | 0.85 (0.75–0.95) |

What I take from this so far:

- Expression carries most of the signal. Mutations alone do much better than chance (1 in 7, about 0.14), but adding them to expression did not improve the scores.
- The simpler logistic regression did at least as well as XGBoost on this small dataset.
- Liver biopsies are a big part of the expression variation (the first principal component largely separates them), but the model was not worse on liver biopsies (out-of-fold accuracy 0.93 vs 0.88 for other sites), and liver-specific genes contributed very little to its predictions.
- Known tissue markers come up among the most important genes, for example GATA3 for breast, CDX2 for colorectal and NKX2-1 for lung.

## Progress

- [x] Repository and test skeleton
- [x] Data loading and cleaning
- [x] Exploratory analysis and PCA
- [x] Baseline model (logistic regression)
- [x] XGBoost and feature-set comparison
- [x] SHAP analysis and biopsy-site check
- [ ] Nextflow pipeline and Docker image tested end to end
- [ ] CI running on GitHub
- [ ] Final report

## Add-on: variant annotation of misclassified cases

**Progress:** worksheets generated; manual curation in progress.

A model can say a tumour looks like pancreatic cancer, but a genetics lab would also check the patient's variants. This add-on asks: **for the patients the model got wrong, do their somatic variants support the true cancer type, the predicted one, or neither?**

`src/annotate_cases.py` takes the 7 misclassified test patients, extracts their protein-changing variants from the raw mutation file (same QC as `clean.py`), and builds an annotation worksheet for the variants in driver genes. `src/gather_evidence.py` then fills in the factual evidence from open databases (CIViC, ClinVar, cBioPortal/TCGA and UniProt). Interpretation is done by hand: each variant is reviewed (with COSMIC checked manually) and given an AMP/ASCO/CAP tier. Each case gets a short review following a lab-style report checklist.

```bash
python src/annotate_cases.py
python src/gather_evidence.py
```

The procedure, evidence sources and templates are in [docs/annotation_guide.md](docs/annotation_guide.md) and [docs/templates/](docs/templates/). Patient-level worksheets and case reviews stay in `results/annotation/`, which is not committed, because the POG570 data can't be redistributed.

First look, before curation: 360 protein-changing variants across the 7 patients, 11 of them in driver genes. Two errors involve tumours whose label groups different histologies, two patients carry variants typical of their true cancer type, and three have rare tumour types, two with no driver-gene variants.

## Limitations

These will be revisited once results are in, but some are clear from the start:

- 570 patients spread over many cancer types means some classes are small. Rare types will be dropped or grouped.
- All samples come from one program at one centre, so results may not hold for other hospitals or sequencing methods.
- Many POG patients had been treated before the biopsy, and the POG570 paper showed that treatment leaves marks on the tumour genome (Pleasance et al., 2020).
- A good score on known cancer types does not mean the model would work on true CUP cases.

## References

- Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*. 2016:785–794. doi:10.1145/2939672.2939785
- Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017;35:316–319. doi:10.1038/nbt.3820
- Grewal JK, Tessier-Cloutier B, Jones M, et al. Application of a neural network whole transcriptome–based pan-cancer method for diagnosis of primary and metastatic cancers. *JAMA Netw Open*. 2019;2(4):e192597. doi:10.1001/jamanetworkopen.2019.2597
- Jiao W, Atwal G, Polak P, et al. A deep learning system accurately classifies primary and metastatic cancers using passenger mutation patterns. *Nat Commun*. 2020;11:728. doi:10.1038/s41467-019-13825-8
- Lundberg SM, Lee SI. A unified approach to interpreting model predictions. *Advances in Neural Information Processing Systems 30 (NIPS 2017)*. 2017:4765–4774.
- Moon I, LoPiccolo J, Baca SC, et al. Machine learning for genetics-based classification and treatment response prediction in cancer of unknown primary. *Nat Med*. 2023;29:2057–2067. doi:10.1038/s41591-023-02482-6
- Pavlidis N, Pentheroudakis G. Cancer of unknown primary site. *Lancet*. 2012;379(9824):1428–1435. doi:10.1016/S0140-6736(11)61178-1
- Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. *J Mach Learn Res*. 2011;12:2825–2830.
- Penson A, Camacho N, Zheng Y, et al. Development of genome-derived tumor type prediction to inform clinical cancer care. *JAMA Oncol*. 2020;6(1):84–91. doi:10.1001/jamaoncol.2019.3985
- Pleasance E, Titmuss E, Williamson L, et al. Pan-cancer analysis of advanced patient tumors reveals interactions between therapy and genomic landscapes. *Nat Cancer*. 2020;1:452–468. doi:10.1038/s43018-020-0050-6
- Sondka Z, Bamford S, Cole CG, et al. The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. *Nat Rev Cancer*. 2018;18:696–705. doi:10.1038/s41568-018-0060-1

## Acknowledgements

This work uses data from the Personalized OncoGenomics (POG) program at BC Cancer. I thank the patients and families who took part, and the POG team for making the data available. The full acknowledgement text requested by the data providers is in `POG570_README.txt`, which comes with the download.

## Author

Maryam Allahyari
