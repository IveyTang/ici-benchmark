# data/

Raw public data live here. **This directory is not tracked by git** — only this
README is. Run `bash scripts/download_data.sh` to populate it.

Expected layout after a successful download:

```
data/
├── xena/                                     TCGA-HNSC (UCSC Xena / GDC hub)
│   ├── TCGA-HNSC.star_fpkm.tsv.gz                ~140 MB
│   ├── TCGA-HNSC.survival.tsv.gz
│   ├── TCGA-HNSC.clinical.tsv.gz
│   └── TCGA-HNSC.somaticmutation_wxs.tsv.gz
├── geo/                                      five ICI validation cohorts
│   ├── GSE159067_series_matrix.txt.gz
│   ├── GSE159067_IHN_log2cpm_data.txt.gz
│   ├── GSE176307_series_matrix.txt.gz
│   ├── GSE176307_norm.csv.gz
│   ├── GSE91061_series_matrix.txt.gz
│   ├── GSE91061_fpkm.csv.gz
│   ├── GSE78220_series_matrix.txt.gz
│   ├── GSE78220_PatientFPKM.xlsx
│   ├── GSE135222_series_matrix.txt.gz
│   └── GSE135222_exp.tsv.gz
├── gh/IMvigor210/data/                       IMvigor210CoreBiologies 1.0.0
│   ├── cds.RData                                 ~51 MB, DESeq v1 CountDataSet
│   ├── dat19.RData  dat25.RData  dat57.RData
│   ├── fmone.RData   color_palettes.RData   datalist
├── msigdb/
│   ├── h.all.v2024.1.Hs.symbols.gmt
│   └── c7.all.v2024.1.Hs.symbols.gmt
└── anno/
    └── Homo_sapiens.gene_info.gz
```

Total: roughly 420 MB.

Per-file provenance, accessions and download URLs:
**[../docs/DATA_SOURCES.md](../docs/DATA_SOURCES.md)**.

Notes
-----
* `cds.RData` holds a **DESeq v1 `CountDataSet`**, not a DESeq2 object. `01_prep.R`
  reads its S4 slots directly, so the retired `DESeq` package is not required —
  do not "fix" this by installing DESeq.
* If NCBI is unreachable from your network, `scripts/ncbi_geo.sh` retries over
  FTP and then HTTPS with an explicit IP resolution. See the network note in
  `docs/DATA_SOURCES.md`.
* `MSigDB` filenames are version-pinned to release 2024.1.Hs. Using a different
  release will change the 110-gene-set panel and therefore the discovery-stage
  phenotype.
* Downloaded files are skipped if already present. Use
  `bash scripts/download_data.sh --force` to refetch.
