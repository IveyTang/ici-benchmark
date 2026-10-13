## ============================================================
## exp_ias.R — IAS 构造方案对照实验
##   比较：全局免疫分数 / 当前 LASSO / 分位 LASSO / PCA1 / 单变量加权
##   评价：6 个 ICI 队列的逐队列 AUC + 随机效应合并 AUC
## ============================================================
suppressPackageStartupMessages({ library(data.table); library(pROC); library(GSVA); library(metafor) })
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))

## ---------- 重建 110 个免疫签名（与 02_discovery.R 完全一致）----------
read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  parts <- strsplit(ln, "\t", fixed = TRUE)
  gs <- lapply(parts, function(x) unique(x[-c(1, 2)]))
  names(gs) <- vapply(parts, function(x) x[1], ""); gs[lengths(gs) > 0]
}
tcX <- readRDS(file.path(RDIR, "tcga_hnsc.rds"))$X
gsc_h <- read_gmt(file.path(DL, "msigdb/h.all.v2024.1.Hs.symbols.gmt"))
c7_all <- read_gmt(file.path(DL, "msigdb/c7.all.v2024.1.Hs.symbols.gmt"))
kw <- paste(c("CD8","CD4","TREG","REGULATORY.T","NK.CELL","NKT","BCELL","B.CELL","DENDRITIC",
              "MACROPHAGE","MONOCYTE","NEUTROPHIL","MAST.CELL","EOSINOPHIL","IFNG","INTERFERON",
              "CYTOLYTIC","TCELL","T.CELL","IMMUNE","GRANZYME","GZMA","GZMB","PD1","PDL1",
              "CHECKPOINT","TH1","TH2","TH17","CENTRAL.MEMORY","EFFECTOR.MEMORY","PLASMA.CELL",
              "MYELOID"), collapse = "|")
gsc_c7 <- c7_all[grepl(kw, names(c7_all), ignore.case = TRUE)]
gsc_c7 <- gsc_c7[lengths(gsc_c7) >= 10 & lengths(gsc_c7) <= 500]
gsc_c7 <- gsc_c7[order(-lengths(gsc_c7))][1:60]
gsc <- c(gsc_h, gsc_c7); gsc <- lapply(gsc, unique); gsc <- gsc[lengths(gsc) >= 10]
gsc <- lapply(gsc, function(g) intersect(g, rownames(tcX))); gsc <- gsc[lengths(gsc) >= 10]
cat("签名数:", length(gsc), "\n")

## ---------- TCGA 与各队列的 GSVA 分数 ----------
TC <- readRDS(file.path(RDIR, "tcga_enrich_gsva.rds"))
coh <- readRDS(file.path(RDIR, "ici_cohorts.rds"))
f_co <- file.path(RDIR, "cohort_enrich_gsva.rds")
if (file.exists(f_co)) { CO <- readRDS(f_co) } else {
  CO <- list()
  for (nm in names(coh)) {
    o <- coh[[nm]]
    gs <- lapply(gsc, function(g) intersect(g, rownames(o$X)))
    ok <- lengths(gs) >= 10
    cat("GSVA:", nm, " 保留签名", sum(ok), "/", length(gsc), "\n")
    CO[[nm]] <- run_enrich(o$X, gs[ok], method = "gsva")
  }
  saveRDS(CO, f_co)
}
for (nm in names(CO)) cat(nm, ":", nrow(CO[[nm]]), "sets x", ncol(CO[[nm]]), "samples\n")

## ---------- 工具函数 ----------
zrow <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }
pooled <- function(auc, n1, n0) {
  d <- data.table(a = auc, n1 = n1, n0 = n0)
  d[, logit := log(a / (1 - a))]
  hm <- function(a, n1, n0) { q1 <- a/(2-a); q2 <- 2*a^2/(1+a)
    (a*(1-a) + (n1-1)*(q1-a^2) + (n0-1)*(q2-a^2))/(n1*n0) }
  d[, vi := hm(a, n1, n0)]
  m <- rma(yi = logit, vi = vi, data = d, method = "DL")
  c(auc = plogis(m$b[1,1]), lo = plogis(m$ci.lb), hi = plogis(m$ci.ub), I2 = m$I2)
}
detail <- function(b, label) {
  zt <- zrow(TC); nb <- intersect(names(b), rownames(zt)); bb <- b[nb] / sum(abs(b[nb]))
  rbindlist(lapply(names(CO), function(nm) {
    zc <- zrow(CO[[nm]]); nn <- intersect(nb, rownames(zc))
    if (length(nn) < 3) return(NULL)
    b2 <- bb[nn] / sum(abs(bb[nn]))
    s <- as.numeric(t(zc[nn, , drop = FALSE]) %*% b2)
    y <- coh[[nm]]$response; keep <- !is.na(y) & is.finite(s)
    data.table(variant = label, cohort = nm,
               auc = as.numeric(pROC::auc(pROC::roc(y[keep], s[keep], quiet = TRUE))))
  }))
}
eval_weights <- function(b, label) {
  dd <- detail(b, label)
  rbindlist(lapply(names(CO), function(nm) {
    zt <- zrow(TC); nb <- intersect(names(b), rownames(zt)); bb <- b[nb] / sum(abs(b[nb]))
    zc <- zrow(CO[[nm]]); nn <- intersect(nb, rownames(zc))
    if (length(nn) < 3) return(NULL)
    b2 <- bb[nn] / sum(abs(bb[nn]))
    s <- as.numeric(t(zc[nn, , drop = FALSE]) %*% b2)
    y <- coh[[nm]]$response; keep <- !is.na(y) & is.finite(s)
    if (sum(keep) < 10 || length(unique(y[keep])) < 2) return(NULL)
    data.table(cohort = nm, auc = as.numeric(pROC::auc(pROC::roc(y[keep], s[keep], quiet = TRUE))),
               n1 = sum(y[keep] == 1), n0 = sum(y[keep] == 0))
  })) -> rows
  p <- pooled(rows$auc, rows$n1, rows$n0)
  data.table(variant = label, n_sig = length(b), mean_auc = mean(rows$auc),
             pooled_auc = p[1], lo = p[2], hi = p[3], I2 = p[4])
}

## ---------- 免疫分数与表型 ----------
zt <- zrow(TC); imm <- colMeans(zt)
lab_lasso <- readRDS(file.path(RDIR, "tcga_phenotype.rds"))$pheno
cat("\n全局免疫分数 sd =", round(sd(imm), 3),
    "; k2 簇间均值差 =", round(diff(range(tapply(imm, lab_lasso, mean))), 3), " SD\n")

v1 <- setNames(rep(1, nrow(zt)), rownames(zt))
v2 <- readRDS(file.path(RDIR, "discovery_results.rds"))$beta
suppressPackageStartupMessages(library(glmnet))
set.seed(20261012)
y3 <- rep(NA_integer_, length(imm)); names(y3) <- names(imm)
y3[imm > quantile(imm, 0.75)] <- 1L; y3[imm < quantile(imm, 0.25)] <- 0L
ok3 <- !is.na(y3)
Xf <- t(zt)                                   # samples x signatures
cv3 <- cv.glmnet(Xf[ok3, ], y3[ok3], family = "binomial", alpha = 1, nfolds = 10)
c3 <- as.matrix(coef(cv3, s = "lambda.1se")); v3 <- c3[c3[,1] != 0 & rownames(c3) != "(Intercept)", 1]
c3b <- as.matrix(coef(cv3, s = "lambda.min")); v3b <- c3b[c3b[,1] != 0 & rownames(c3b) != "(Intercept)", 1]
pc <- prcomp(Xf, center = FALSE, scale. = FALSE)
v4 <- pc$rotation[, 1]; names(v4) <- rownames(zt)
if (cor(as.numeric(Xf %*% v4), imm) < 0) v4 <- -v4
tt <- apply(Xf[ok3, ], 2, function(x) t.test(x[y3[ok3] == 1], x[y3[ok3] == 0])$statistic)
v5 <- setNames(sort(tt, decreasing = TRUE)[1:20], names(sort(tt, decreasing = TRUE)[1:20]))
cr <- as.numeric(cor(Xf, imm)); names(cr) <- rownames(zt)
v6 <- setNames(sort(cr, decreasing = TRUE)[1:20], names(sort(cr, decreasing = TRUE)[1:20]))

cat("V3(1se) 特征数:", length(v3), " V3b(min) 特征数:", length(v3b), "\n")
res <- rbindlist(list(
  eval_weights(v1, "V1_global_immune"),
  eval_weights(v2, "V2_lasso_k2"),
  eval_weights(v3, "V3_lasso_quartile_1se"),
  eval_weights(v3b, "V3b_lasso_quartile_min"),
  eval_weights(v4, "V4_pca1"),
  eval_weights(v5, "V5_univ_top20"),
  eval_weights(v6, "V6_corr_top20")
))
cat("\n===== 方案对照 =====\n"); print(res, digits = 3)
cat("\n===== 逐队列 AUC =====\n")
dd <- rbindlist(list(list(v1,"V1_global_immune"), list(v2,"V2_lasso_k2"),
  list(v3,"V3_lasso_quartile_1se"), list(v3b,"V3b_lasso_quartile_min"), list(v4,"V4_pca1"),
  list(v5,"V5_univ_top20"), list(v6,"V6_corr_top20")) |> lapply(function(z) detail(z[[1]], z[[2]])))
print(dcast(dd, variant ~ cohort, value.var = "auc"), digits = 3)
saveRDS(list(res = res, detail = dd, v3 = v3, v3b = v3b), file.path(RDIR, "exp_ias_results.rds"))
