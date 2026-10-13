## ============================================================
## 07_figures_extra.R — 基准比较森林图（fig6）+ 免疫细胞标记基因面板（fig7）
## ============================================================
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))

e1 <- readRDS(file.path(RDIR, "exp_ias_results.rds"))$res
e2 <- readRDS(file.path(RDIR, "exp_ias2_results.rds"))$R
e3 <- readRDS(file.path(RDIR, "exp_ias3_results.rds"))$R
e4 <- readRDS(file.path(RDIR, "exp_ias4_results.rds"))$S

lab <- c(
  V1_global_immune         = "Global immune score (GSVA, 110 sets)",
  V2_lasso_k2              = "LASSO on GSVA scores (40 sets)",
  V3_lasso_quartile_1se    = "LASSO on GSVA scores, quartile (18 sets)",
  V3b_lasso_quartile_min   = "LASSO on GSVA scores, quartile (20 sets)",
  V4_pca1                  = "PC1 of GSVA scores (110 sets)",
  V5_univ_top20            = "Univariate top-20 GSVA sets",
  V6_corr_top20            = "Correlation top-20 GSVA sets",
  L_pheno_1se              = "LASSO, expression, phenotype label (69 genes)",
  L_pheno_min              = "LASSO, expression, phenotype label (110 genes)",
  L_quart_1se              = "LASSO, expression, immune-score quartile (48 genes)",
  L_quart_min              = "LASSO, expression, immune-score quartile (51 genes)",
  U_top20                  = "Univariate top-20 genes",
  U_top50                  = "Univariate top-50 genes",
  IMMUNE_ALL               = "All immune genes (9,661)",
  COR_top10                = "Correlation-weighted top-10 genes",
  COR_top15                = "Correlation-weighted top-15 genes",
  COR_top20                = "Correlation-weighted top-20 genes",
  COR_top30                = "Correlation-weighted top-30 genes",
  EN_a0.5_1se              = "Elastic net (a=0.5, IFN-g axis, 109 genes)",
  EN_a0.5_min              = "Elastic net (a=0.5, IFN-g axis, 110 genes)",
  EN_a1.0_1se              = "LASSO (IFN-g axis, 26 genes)  =  IAS",
  EN_a1.0_min              = "LASSO (IFN-g axis, 27 genes)",
  REF_IFNG6                = "IFN-g 6-gene signature (Ayers)",
  REF_AYERS                = "Ayers 19-gene signature",
  REF_IFNG_HALL            = "Hallmark IFN-g response (199 genes)",
  REF_CYT2                 = "Cytolytic activity (GZMA, PRF1)",
  F_IFNG_6                 = "IFN-g 6-gene signature",
  F_TIS_18                 = "TIS 18-gene signature",
  TRAIN_IMV                = "LASSO trained in IMvigor210 (8 genes)")

A <- rbindlist(list(
  e1[, .(variant, pooled_auc, lo, hi)],
  e2[, .(variant, pooled_auc, lo, hi)],
  e3[, .(variant, pooled_auc, lo, hi)],
  e4[, .(variant = ifelse(scheme == "TrainIMvigor210", "TRAIN_IMV", scheme),
         pooled_auc, lo, hi)][grepl("^REF_|^TRAIN_", variant)]))
A <- unique(A, by = "variant")
A <- A[!variant %in% c("F_IFNG_6", "REF_CYT2")]
A[, kind := fifelse(grepl("^REF_|^TRAIN_", variant), "Published / reference",
             fifelse(grepl("^F_", variant), "Published / reference", "Derived in this study"))]
A[, group := fifelse(grepl("^V[0-9]", variant), "Gene-set enrichment score level",
             fifelse(grepl("^REF_|^TRAIN_|^F_", variant), "Reference signatures",
                     "Gene-expression level"))]
setorder(A, pooled_auc)
A[, label := ifelse(is.na(lab[variant]), variant, lab[variant])]

png(file.path(FDIR, "fig6_benchmark_strategies.png"), 3600, 2300, res = 300)
par(mar = c(4.5, 26, 3.5, 5))
colg <- c("Gene-set enrichment score level" = "#2c5aa0",
          "Gene-expression level" = "#c0392b",
          "Reference signatures" = "#27ae60")
plot(A$pooled_auc, seq_len(nrow(A)), xlim = c(0.48, 0.72), pch = 18, cex = 1.6,
     col = colg[A$group], yaxt = "n", xlab = "Pooled AUC (random-effects, 6 cohorts)",
     ylab = "", main = "Cross-cohort discrimination by score-construction strategy")
segments(A$lo, seq_len(nrow(A)), A$hi, seq_len(nrow(A)), lwd = 1.8, col = colg[A$group])
abline(v = 0.5, lty = 3, col = "grey50")
abline(v = 0.6676, lty = 2, col = "#27ae60", lwd = 2)
axis(2, at = seq_len(nrow(A)), labels = A$label, las = 1, cex.axis = 0.72, tick = FALSE)
legend("bottomright", names(colg), col = colg, lwd = 2.5, bty = "n", cex = 0.85)
dev.off()
msg("fig6_benchmark_strategies.png 完成")

## ---- fig7: 免疫细胞标记基因面板（替代 LM22 反卷积）----
f_cmp <- file.path(TDIR, "marker_vs_phenotype.csv")
f_ias <- file.path(TDIR, "marker_vs_ias.csv")
if (file.exists(f_cmp)) {
  cmp <- fread(f_cmp)
  setorder(cmp, cohens_d)
  cmp[, star := fifelse(fdr < 0.001, "***", fifelse(fdr < 0.01, "**",
                    fifelse(fdr < 0.05, "*", "ns")))]

  rho <- if (file.exists(f_ias)) fread(f_ias) else NULL

  png(file.path(FDIR, "fig7_cell_marker_panel.png"), 3700, 1900, res = 300)
  par(mfrow = c(1, 2), mar = c(4.5, 12, 3.8, 3), mgp = c(2.6, 0.8, 0))

  ## 面板 A：Cohen's d
  colA <- ifelse(cmp$cohens_d > 0, "#c0392b", "#2c5aa0")
  bp <- barplot(cmp$cohens_d, horiz = TRUE, col = colA, border = NA,
                names.arg = cmp$cell_type, las = 1, cex.names = 0.85,
                xlab = "Cohen's d  (activated vs cold)",
                main = "A  Marker scores by immune phenotype", cex.main = 1.0,
                xlim = c(min(0, min(cmp$cohens_d)) - 0.12, max(cmp$cohens_d) + 0.22))
  abline(v = 0, col = "grey40")
  text(cmp$cohens_d + ifelse(cmp$cohens_d > 0, 0.03, -0.03), bp,
       cmp$star, cex = 0.9, col = "grey20",
       pos = ifelse(cmp$cohens_d > 0, 4, 2))

  ## 面板 B：与 IAS 的 Spearman 相关
  if (!is.null(rho) && nrow(rho)) {
    setorder(rho, rho_IAS)
    colB <- ifelse(rho$rho_IAS > 0, "#c0392b", "#2c5aa0")
    barplot(rho$rho_IAS, horiz = TRUE, col = colB, border = NA,
            names.arg = rho$cell_type, las = 1, cex.names = 0.85,
            xlab = expression("Spearman " * rho * " with IAS"),
            main = "B  Concordance with the derived score", cex.main = 1.0,
            xlim = c(0, max(rho$rho_IAS) * 1.18))
    abline(v = 0, col = "grey40")
    text(rho$rho_IAS + max(rho$rho_IAS) * 0.02, seq_len(nrow(rho)),
         sprintf("%.2f", rho$rho_IAS), cex = 0.78, pos = 4, col = "grey20")
  } else {
    plot.new(); title(main = "B  (IAS concordance unavailable)")
  }
  dev.off()
  msg("fig7_cell_marker_panel.png 完成")
  print(cmp)
} else msg("标记面板结果尚未生成，跳过 fig7")
