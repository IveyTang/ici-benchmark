# Data sources

All data analysed in this study are publicly available. No new primary data were
generated. Raw data are **not** redistributed in this repository; the download
script `scripts/download_data.sh` retrieves every input file into `data/`, and
`scripts/01_prep.R` rebuilds all derived objects from it.

---

## 1. Discovery cohort — TCGA-HNSC

Source: UCSC Xena / GDC harmonised hub (`gdc-hub.s3.us-east-1.amazonaws.com`).

| File placed in `data/xena/` | Content |
|---|---|
| `TCGA-HNSC.star_fpkm.tsv.gz` | STAR-FPKM expression, 60,483 features x 566 samples |
| `TCGA-HNSC.survival.tsv.gz` | OS time (days) and OS event |
| `TCGA-HNSC.clinical.tsv.gz` | GDC clinical attributes |
| `TCGA-HNSC.somaticmutation_wxs.tsv.gz` | WXS somatic mutations (used for TMB) |

Only primary tumours were retained (sample barcode positions 14-15 equal to
`01`), giving **520 tumours** and **221 death events**. TMB was computed as
`n_mutations / 38 Mb` and is therefore a study-specific approximation, not the
TCGA-published call set.

Also used: `data/xena/TCGA-HNSC.GDC_phenotype.tsv.gz` is not required by any
script and is not downloaded.

## 2. ICI-treated validation cohorts

| Cohort | Cancer / therapy | Samples analysed (responders) | Accession |
|---|---|---|---|
| GSE159067 | HNSCC, anti-PD-1 | 102 (11) | GEO |
| GSE176307 | Urothelial (BACI), anti-PD-1 | 88 (16) | GEO |
| IMvigor210 | Urothelial, atezolizumab | 298 assessable (68) | see below |
| GSE91061 | Melanoma, nivolumab | 49 (10) | GEO |
| GSE78220 | Melanoma, pembrolizumab | 27 (15) | GEO |
| GSE135222 | NSCLC, anti-PD-1/PD-L1 | 27 (7) | GEO |

Total **591 patients, 127 responders**.

### 2.1 GEO cohorts (GSE159067, GSE176307, GSE91061, GSE78220, GSE135222)

For each series the analysis requires the **series matrix** (sample annotation)
plus the **supplementary expression file**:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSExxxnnn/GSExxxxx/matrix/GSExxxxx_series_matrix.txt.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSExxxnnn/GSExxxxx/suppl/<expression file>
```

where `GSExxxnnn` is the accession with the last three digits replaced by
`nnn` (e.g. GSE159067 -> `GSE159nnn`). `scripts/ncbi_geo.sh` resolves these
paths and downloads them; see "Network note" below.

Expression files used:

| Accession | Supplementary file | Platform |
|---|---|---|
| GSE159067 | `GSE159067_IHN_log2cpm_data.txt.gz` | HTG EdgeSeq (log2 CPM) |
| GSE176307 | `GSE176307_norm.csv.gz` | RNA-seq (normalised) |
| GSE91061 | `GSE91061_fpkm.csv.gz` | RNA-seq (FPKM) |
| GSE78220 | `GSE78220_PatientFPKM.xlsx` | RNA-seq (FPKM) |
| GSE135222 | `GSE135222_exp.tsv.gz` | RNA-seq |

Sample annotation is parsed cell-by-cell rather than row-wise, because GEO
`!Sample_characteristics` rows are misaligned when samples carry different
numbers of characteristics (see the comment on `series_matrix_meta()` in
`01_prep.R`). For GSE91061 only baseline (pre-treatment) samples were used.

### 2.2 IMvigor210

The IMvigor210 cohort (Mariathasan et al., *Nature* 2018) is distributed as the
R data package **IMvigor210CoreBiologies** under CC-BY 3.0:

```
http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz
```

The underlying sequencing data are deposited at the European Genome-Phenome
Archive under accession **EGAS00001002556**.

Important implementation detail: `data/cds.RData` contains a **DESeq v1
`CountDataSet`**, not a DESeq2 or SummarizedExperiment object. Methods such as
`assay()`, `exprs()` and `as(cds, "ExpressionSet")` therefore fail unless the
retired `DESeq` package is installed. `01_prep.R` reads the S4 slots directly:

```r
M <- get("counts", envir = attr(cds, "assayData"))        # Entrez x sample
pdi <- as.data.frame(attr(attr(cds, "phenoData"), "data"))
```

Response was taken from `binaryResponse` (`CR/PR` -> responder, `SD/PD` ->
non-responder). Overall survival is stored in `os` **in months** and is
converted to days (`x 30.4375`) for consistency with the GEO cohorts; `censOS`
equal to 1 marks a death event.

## 3. Gene sets and annotation

| File placed in `data/` | Source |
|---|---|
| `msigdb/h.all.v2024.1.Hs.symbols.gmt` | MSigDB Hallmark, release 2024.1.Hs |
| `msigdb/c7.all.v2024.1.Hs.symbols.gmt` | MSigDB C7 immunologic signatures |
| `anno/Homo_sapiens.gene_info.gz` | NCBI Gene `GENE_INFO/Mammalia/` |

MSigDB files come from `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/`.
The manuscript states that **110 immune gene sets** were used: this is the C7
collection after restricting to sets with at least 10 and at most 500 genes
present in the TCGA-HNSC matrix, and after removing sets overlapping the
IFN-gamma axis used to define the anchor label.

The 6-gene IFN-gamma signature (IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG), the
Ayers 19-gene signature and the TIS 18-gene signature are hard-coded in
`06_final.R`; the Danaher marker panel is hard-coded in `04_markers.R`.

---

## Network note (NCBI / GEO accessibility)

On some networks `ftp.ncbi.nlm.nih.gov` and `www.ncbi.nlm.nih.gov` are blocked
at the domain level, while FTP and direct-IP HTTPS remain reachable.
`scripts/ncbi_geo.sh` therefore tries two channels in order:

1. **FTP** — `ftp://ftp.ncbi.nlm.nih.gov/geo/series/...` (preferred, simplest).
2. **HTTPS with explicit resolution** —
   `curl --resolve ftp.ncbi.nlm.nih.gov:443:<IP> --noproxy '*' -k https://...`

If your network allows ordinary HTTPS to NCBI, plain `curl`/`wget` will also
work and the script's fallbacks are harmless. If NCBI is entirely unreachable,
the same series matrices are mirrored by EBI ArrayExpress and the expression
matrices can be obtained from the corresponding publications' supplementary
files.
