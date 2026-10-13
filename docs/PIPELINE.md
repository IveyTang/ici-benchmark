# Pipeline

Scripts are numbered in execution order. All of them resolve their paths
through `scripts/00_config.R`, so they can be run from the repository root
without editing anything:

```bash
export ICI_WORK="$(pwd)"      # already the default when you run from the root
Rscript scripts/01_prep.R
```

`bash scripts/run_all.sh` runs the whole chain and stops at the first failure.

---

## Stage 0 — data acquisition

`scripts/download_data.sh` (and `scripts/ncbi_geo.sh` for the GEO series).
See `docs/DATA_SOURCES.md`. Output goes to `data/`.

## Stage 1 — `01_prep.R`: cohorts and expression matrices

* Builds the Entrez/Ensembl to symbol maps from `Homo_sapiens.gene_info.gz`,
  collapsing duplicate symbols by the mean of their rows.
* **TCGA-HNSC**: keeps primary tumours only (`01` in barcode positions 14-15),
  drops features with any non-finite value or zero variance; normalises the
  counts-free FPKM matrix to log2(CPM + 1); derives TMB from WXS mutations
  (`n / 38 Mb`).
* **Five GEO cohorts** and **IMvigor210** are harmonised onto a common symbol
  space and stored as a list of `list(X, meta)` objects.
* Parses GEO series matrices cell-by-cell. This matters: `!Sample_characteristics`
  rows are misaligned whenever samples carry different numbers of
  characteristics, so a row-wise split silently shifts values between fields.
* All heavy steps are cached under `results/rds/`; delete the cache to force a
  rebuild.

Expected inclusions (verified against the frozen run): GSE159067 102/11,
GSE78220 27/15, GSE91061 49/10, GSE135222 27/7, GSE176307 88/16,
IMvigor210 348 samples of which 298 assessable (68/230).

## Stage 2 — `02_discovery.R`: immune phenotype in TCGA-HNSC

* GSVA enrichment scores (`kcdf = "Gaussian"`, `maxDiff = TRUE`) over the 110
  immune gene sets → `results/rds/tcga_enrich_gsva.rds`.
* Consensus clustering of the **samples** over those scores, `maxK = 6`,
  `reps = 1000`, Spearman distance, average linkage. Output `fig0`,
  `tcga_phenotype.rds`.
* Derivative score (`IAS`): a 700-gene pool is defined around the IFN-gamma
  axis, and `glmnet` logistic LASSO (`family = "binomial"`, `lambda.1se`) is
  fitted with the IFN-gamma-axis anchor label, yielding a **26-gene** model.
  Coefficient table → `results/tables/ias_coefficients.csv`.

> **Why GSVA and not ssGSEA here.** Both are within-sample rank-based
> statistics. With ssGSEA the median between-sample Spearman correlation of the
> enrichment matrix is approximately **0.93**, and consensus clustering
> degenerates every time to 517/3. With GSVA the same correlation is
> approximately **0**, and k = 2 returns a balanced 251/269 split. The
> degeneracy is a property of the statistic, not of the data.

> **Why the anchor is independent of the label.** If the phenotype label were
> the median split of the features themselves, the training AUC would be 1.0 by
> construction. The label is therefore derived from an independent IFN-gamma
> axis, and a 70/30 stratified split separates fitting from evaluation.

## Stage 3 — exploration: `exp_ias.R`, `exp_ias2.R`, `exp_ias3.R`, `exp_ias4.R`

Four independent experiments that together define the **27 score-construction
strategies** of the benchmark. They are kept as separate scripts because each
answers a different question, and because collapsing them would hide which
strategy contributed which variant.

| Script | Question | Result object |
|---|---|---|
| `exp_ias.R` | Set-level aggregation: how should gene-set scores be combined into one number per patient? (7 variants) | `exp_ias_results.rds` |
| `exp_ias2.R` | Gene-level aggregation with different label definitions (9 variants) | `exp_ias2_results.rds` |
| `exp_ias3.R` | Variants derived from the IFN-gamma anchor (11 variants) | `exp_ias3_results.rds` |
| `exp_ias4.R` | Cross-cohort leave-one-out training, including training inside a validation cohort (5 schemes) | `exp_ias4_results.rds` |

All four use the identical scoring, discrimination and pooling code so that
their outputs are directly comparable.

## Stage 4 — `04_markers.R`: immune cell marker panel

Danaher et al. marker genes for 13 lineages; per-gene z-scoring across samples,
then the mean per lineage. Reports Cohen's d (Hedges-corrected) between immune
phenotypes and the Spearman correlation of each lineage with the derived score.

> **Rejected approach.** Signature-matrix deconvolution with the LM22 matrix was
> attempted first and abandoned. Under the published CIBERSORT
> parameterisation (nu-support vector regression, `scale = FALSE`) a single
> sample failed to converge in over three minutes, and six parallel workers
> produced no result in twenty minutes. With `scale = TRUE`, or with the
> signature columns manually standardised, fitting became fast but the median
> correlation between fitted and observed profiles was only 0.05, against the
> 0.4-0.9 expected of a well-conditioned deconvolution. Constrained least
> squares (`quadprog`, equality plus non-negativity) reported inconsistent
> constraints; `nnls::nnls` reached a correlation of 0.60 but returned
> degenerate solutions in which two of twenty-two lineages absorbed almost all
> of the weight (`B cells memory` alone at 0.814). The cause is intrinsic: LM22
> is highly co-linear and is not quantile-normalised, so a convex combination
> cannot approximate bulk head and neck tumour profiles. The marker panel is
> used instead, and it is also more consistent with the paper's argument that
> transparent expression-level aggregation is preferable.

## Stage 5 — `06_final.R`: cross-cohort validation (the primary results)

* Scores every strategy in all six ICI cohorts, z-scored **within cohort**.
* Per-cohort AUC with DeLong confidence intervals; odds ratio per standard
  deviation; survival hazard ratio where follow-up is available.
* Pooled AUC by **DerSimonian-Laird random effects on the logit scale**, with
  Hanley-McNeil variance and a floor/cap on the logit to keep degenerate
  cohorts finite.
* Leave-one-cohort-out (LOCO) and a permutation test (200 permutations of the
  gene labels) for the derived score.
* Produces `fig1`-`fig5`, `final_results.rds`, and tables
  `final_cohort_performance.csv`, `final_pooled.csv`,
  `final_leave_one_out.csv`, `final_ias_scores.csv`,
  `table1_tcga_baseline.csv`, `table2_cohort_performance.csv`,
  `table3_signature_comparison.csv`.

Headline numbers from the frozen run: IFN-gamma 6-gene 0.668 (95% CI
0.573-0.750); IFN-gamma 6-gene (TIS-style) 0.638; Ayers 19-gene 0.637;
derived 26-gene score 0.601 (95% CI 0.540-0.659; OR/SD 1.37, P = 0.0030).

## Stage 6 — `07_figures_extra.R`: benchmark and marker figures

* `fig6` — forest plot of the 27 strategies, coloured by aggregation level.
* `fig7` — the marker panel (Cohen's d and concordance with the derived score).
  Requires `marker_vs_phenotype.csv` / `marker_vs_ias.csv` from stage 4.

## Stage 7 — `09_benchmark_table.R`: the 27-strategy benchmark table

Rebuilds `results/tables/table4_benchmark_all_strategies.csv` from the four
`exp_ias*_results.rds` objects. This script reproduces the shipped CSV
byte-for-byte, and its level summary is the direct evidence for the paper's
central claim:

| Level | pooled AUC range | n strategies |
|---|---|---|
| Gene expression | 0.5194 - 0.6676 | 19 |
| Published reference / trained | 0.5731 | 1 |
| Gene-set enrichment score | 0.5249 - 0.5571 | 7 |

## Stage 8 — `08_jtm_tables.R`: tables for the manuscript

Extracts the real LASSO coefficients and completes the baseline table with
P values computed from the phenotype assignment (`jtm_table1_baseline_pvalues.csv`,
`jtm_table2_ias_genes.csv`). Verified to reproduce the shipped CSVs
byte-for-byte.

---

## Reproducibility status

`09_benchmark_table.R` and `08_jtm_tables.R` were re-run against the frozen
intermediates from a clean checkout layout and produced output **identical in
content** to the files in `results/tables/` (compared after line-ending
normalisation: the repository stores LF, R's `fwrite` emits CRLF on Windows, and
`.gitattributes` enforces the normalisation). Stages 1-6 involve long-running
steps (the GSVA over 110 sets, the 1000-replicate consensus clustering, and the
permutation test) and were executed once; their outputs are archived in
`results/`.

## Known limitations

* `GSE91061` contributes 49 baseline samples and 10 responders; its per-cohort
  AUC is correspondingly unstable and it is down-weighted by the random-effects
  model but not otherwise corrected.
* TMB is a study-side count rather than a harmonised call set, so its
  association with the derived score should be read as descriptive.
* The 26-gene score is derived in TCGA-HNSC; there is no independent
  derivation cohort, so its coefficients carry no validation of the
  feature-selection step itself, only of the resulting score.
