## ============================================================
## 04_markers.R — 免疫细胞标记基因面板（Danaher 2017 等）
##   替代 LM22 反卷积：透明、确定性、无解退化
##   每个细胞类型 = 其标记基因的 z 化均值表达
## ============================================================
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))

MARK <- list(
  `CD8 T cells`        = c("CD8A", "CD8B"),
  `Cytotoxic cells`    = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7", "KLRD1", "GNLY"),
  `T cells (general)`  = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "CD247"),
  `Regulatory T cells` = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  `B cells`            = c("CD19", "MS4A1", "CD79A", "CD79B", "TNFRSF17", "TCL1A"),
  `NK cells`           = c("KLRF1", "NCR1", "KIR2DL3", "KLRB1", "SH2D1B"),
  `M1 macrophages`     = c("NOS2", "CXCL9", "CXCL10", "IDO1", "IL1B", "SOCS3"),
  `M2 macrophages`     = c("CD163", "MSR1", "MRC1", "CD209", "IL10", "MAF"),
  `Dendritic cells`    = c("CD1C", "CLEC9A", "CCL19", "HSD11B1", "LAMP3"),
  `Mast cells`         = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT"),
  `Neutrophils`        = c("FCGR3B", "CSF3R", "CEACAM8", "S100A12"),
  `Endothelial cells`  = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  `Fibroblasts`        = c("COL1A1", "COL1A2", "DCN", "FAP", "LUM"))

tc    <- readRDS(file.path(RDIR, "tcga_hnsc.rds"))
pheno <- readRDS(file.path(RDIR, "tcga_phenotype.rds"))$pheno
X     <- tc$X                                   # 基因 x 样本，log2(FPKM+1)

miss <- setdiff(unique(unlist(MARK)), rownames(X))
if (length(miss)) msg("警告：以下标记基因不在表达矩阵中 -> ", paste(miss, collapse = ", "))

## 逐基因 z 化（跨样本），再对每个细胞类型的标记取均值
z <- t(scale(t(X)))
score_ct <- function(gs) {
  gs <- intersect(gs, rownames(z))
  if (!length(gs)) return(rep(NA_real_, ncol(z)))
  colMeans(z[gs, , drop = FALSE], na.rm = TRUE)
}
S <- vapply(MARK, score_ct, numeric(ncol(z)))
colnames(S) <- names(MARK)
S <- S[, colSums(is.na(S)) < nrow(S) / 2, drop = FALSE]
msg("标记面板：", ncol(S), " 个细胞类型，", nrow(S), " 个样本")
msg("各细胞类型实际使用基因数：")
for (ct in colnames(S)) msg("  ", sprintf("%-19s", ct),
                        length(intersect(MARK[[ct]], rownames(X))), " / ", length(MARK[[ct]]))

saveRDS(S, file.path(RDIR, "marker_panel.rds"))
fwrite(as.data.table(S, keep.rownames = "sample"), file.path(TDIR, "marker_panel.csv"))

## ---- 与免疫表型比较（Wilcoxon + BH） ----
m  <- data.table(sample = rownames(S), S)
m  <- merge(m, data.table(sample = names(pheno), pheno = as.integer(pheno)), by = "sample")
msg("标记面板-表型合并 n = ", nrow(m))

cmp <- rbindlist(lapply(colnames(S), function(ct) {
  v <- m[[ct]]
  if (all(is.na(v))) return(NULL)
  w <- tryCatch(wilcox.test(v ~ pheno, data = m), error = function(e) NULL)
  if (is.null(w)) return(NULL)
  a <- v[m$pheno == 1]; b <- v[m$pheno == 0]
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  d  <- if (is.finite(sp) && sp > 0) (mean(a) - mean(b)) / sp else NA_real_
  data.table(cell_type = ct,
             med_activated = median(m[pheno == 1][[ct]], na.rm = TRUE),
             med_cold      = median(m[pheno == 0][[ct]], na.rm = TRUE),
             cohens_d      = d, p = w$p.value)
}))
if (!nrow(cmp)) cmp <- data.table(cell_type = character(0), med_activated = numeric(0),
                                  med_cold = numeric(0), cohens_d = numeric(0), p = numeric(0))
cmp[, fdr := p.adjust(p, "BH")]
setorder(cmp, p)
fwrite(cmp, file.path(TDIR, "marker_vs_phenotype.csv"))

## ---- 与既有免疫签名的一致性（验证面板生物学合理性）----
fsc <- file.path(TDIR, "tcga_ias_scores.csv")
if (file.exists(fsc)) {
  fdt <- fread(fsc)                     # sample, IAS, phenotype, split
  fdt <- fdt[!is.na(IAS) & !is.na(phenotype)]
  j <- merge(as.data.table(S, keep.rownames = "sample"),
             fdt[, .(sample, IAS, phenotype)], by = "sample")
  msg("一致性检查可用样本 n = ", nrow(j), "（IAS 与表型标签均非缺失）")
  if (nrow(j) > 20) {
    cr <- rbindlist(lapply(colnames(S), function(ct) {
      ct2 <- tryCatch(cor.test(j[[ct]], j$IAS, method = "spearman"), error = function(e) NULL)
      if (is.null(ct2)) return(NULL)
      data.table(cell_type = ct, rho_IAS = unname(ct2$estimate), p = ct2$p.value)
    }))
    if (nrow(cr)) {
      cr[, fdr := p.adjust(p, "BH")]; setorder(cr, -rho_IAS)
      fwrite(cr, file.path(TDIR, "marker_vs_ias.csv"))
      msg("与本研究 IAS 相关性最高的细胞类型：")
      print(head(cr, 8))
    }
  }
}

msg("标记面板差异最显著的细胞类型：")
print(head(cmp, 8))
msg("MARKER PANEL DONE")
