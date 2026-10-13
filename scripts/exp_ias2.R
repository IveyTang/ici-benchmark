## ============================================================
## exp_ias2.R — 基因表达级 IAS 构造实验
##   结论驱动：GSVA 富集分数丢失免疫丰度信号，改用基因表达 z 分数加权
## ============================================================
suppressPackageStartupMessages({ library(data.table); library(pROC); library(metafor); library(glmnet) })
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))

tc <- readRDS(file.path(RDIR, "tcga_hnsc.rds")); TX <- tc$X
coh <- readRDS(file.path(RDIR, "ici_cohorts.rds"))
zc <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }

## ---------- 参考签名 ----------
sig_fixed <- list(
  IFNG_6 = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"),
  CYT_2  = c("GZMA","PRF1"),
  TIS_18 = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","HLA-DQA1","HLA-DRB1",
             "HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT","TAP1"))

## ---------- TCGA 免疫表型（沿用 GSVA 聚类结果作为「免疫活化」标签）----------
lab <- readRDS(file.path(RDIR, "tcga_phenotype.rds"))$pheno
cat("表型: "); print(table(lab))

## ---------- 免疫基因全集（110 个签名的并集 ∩ TCGA）----------
read_gmt <- function(f) { ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  p <- strsplit(ln, "\t", fixed = TRUE)
  gs <- lapply(p, function(x) unique(x[-c(1, 2)])); names(gs) <- vapply(p, function(x) x[1], ""); gs }
gsc_h <- read_gmt(file.path(DL, "msigdb/h.all.v2024.1.Hs.symbols.gmt"))
c7 <- read_gmt(file.path(DL, "msigdb/c7.all.v2024.1.Hs.symbols.gmt"))
kw <- paste(c("CD8","CD4","TREG","NK.CELL","BCELL","B.CELL","DENDRITIC","MACROPHAGE","MONOCYTE",
  "NEUTROPHIL","MAST.CELL","EOSINOPHIL","IFNG","INTERFERON","CYTOLYTIC","TCELL","T.CELL","IMMUNE",
  "GRANZYME","GZMA","GZMB","PD1","PDL1","CHECKPOINT","TH1","TH2","TH17","PLASMA.CELL","MYELOID"),
  collapse = "|")
g7 <- c7[grepl(kw, names(c7), ignore.case = TRUE)]; g7 <- g7[lengths(g7) >= 10 & lengths(g7) <= 500]
g7 <- g7[order(-lengths(g7))][1:60]
immune_genes <- unique(unlist(c(gsc_h, g7)))
immune_genes <- intersect(immune_genes, rownames(TX))
cat("免疫基因全集:", length(immune_genes), "\n")

## ---------- TCGA 特征矩阵 & 标签 ----------
G <- intersect(unique(c(immune_genes, unlist(sig_fixed))), rownames(TX))
ZT <- zc(TX[G, ])                                   # genes x samples
Xf <- t(ZT)
y <- lab[rownames(Xf)]
imm_exp <- colMeans(ZT[immune_genes, ])             # 表达级全局免疫分数
names(imm_exp) <- colnames(TX)

set.seed(20261012)
yq <- rep(NA_integer_, length(imm_exp)); names(yq) <- names(imm_exp)
yq[imm_exp > quantile(imm_exp, .75)] <- 1L; yq[imm_exp < quantile(imm_exp, .25)] <- 0L

## ---------- 各方案的基因权重 ----------
W <- list()
for (nm in names(sig_fixed)) { g <- intersect(sig_fixed[[nm]], rownames(ZT))
  W[[paste0("F_", nm)]] <- setNames(rep(1, length(g)), g) }
cvA <- cv.glmnet(Xf, y, family = "binomial", alpha = 1, nfolds = 10)
for (s in c("lambda.1se", "lambda.min")) {
  cc <- as.matrix(coef(cvA, s = s)); v <- cc[cc[,1] != 0 & rownames(cc) != "(Intercept)", 1]
  W[[paste0("L_pheno_", sub("lambda\\.", "", s))]] <- v
}
cvB <- cv.glmnet(Xf[!is.na(yq), ], yq[!is.na(yq)], family = "binomial", alpha = 1, nfolds = 10)
for (s in c("lambda.1se", "lambda.min")) {
  cc <- as.matrix(coef(cvB, s = s)); v <- cc[cc[,1] != 0 & rownames(cc) != "(Intercept)", 1]
  W[[paste0("L_quart_", sub("lambda\\.", "", s))]] <- v
}
tt <- apply(Xf, 2, function(x) t.test(x[y == 1], x[y == 0])$statistic)
W[["U_top20"]] <- setNames(sort(tt, decreasing = TRUE)[1:20], names(sort(tt, decreasing = TRUE))[1:20])
W[["U_top50"]] <- setNames(sort(tt, decreasing = TRUE)[1:50], names(sort(tt, decreasing = TRUE))[1:50])
W[["IMMUNE_ALL"]] <- setNames(rep(1, length(immune_genes)), immune_genes)

cat("\n各方案基因数:\n"); print(sapply(W, length))

## ---------- 在 6 个队列上评价 ----------
pooled <- function(a, n1, n0) {
  d <- data.table(a = a, n1 = n1, n0 = n0); d[, logit := log(a/(1-a))]
  hm <- function(a,n1,n0){q1<-a/(2-a);q2<-2*a^2/(1+a);(a*(1-a)+(n1-1)*(q1-a^2)+(n0-1)*(q2-a^2))/(n1*n0)}
  d[, vi := hm(a, n1, n0)]; m <- rma(yi = logit, vi = vi, data = d, method = "DL")
  c(auc = plogis(m$b[1,1]), lo = plogis(m$ci.lb), hi = plogis(m$ci.ub), I2 = m$I2)
}
ZC <- lapply(names(coh), function(nm) zc(coh[[nm]]$X)); names(ZC) <- names(coh)
detail <- function(bv, label) {
  rbindlist(lapply(names(coh), function(nm) {
    g <- intersect(names(bv), rownames(ZC[[nm]])); if (length(g) < 3) return(NULL)
    b <- bv[g] / sum(abs(bv[g]))
    s <- as.numeric(t(ZC[[nm]][g, , drop = FALSE]) %*% b)
    yv <- coh[[nm]]$response; keep <- !is.na(yv) & is.finite(s)
    if (sum(keep) < 10 || length(unique(yv[keep])) < 2) return(NULL)
    data.table(variant = label, cohort = nm, n_gene_used = length(g),
      coverage = length(g)/length(bv),
      auc = as.numeric(pROC::auc(pROC::roc(yv[keep], s[keep], quiet = TRUE))),
      n1 = sum(yv[keep] == 1), n0 = sum(yv[keep] == 0))
  }))
}
D <- rbindlist(lapply(names(W), function(nm) detail(W[[nm]], nm)))
R <- D[, { p <- pooled(auc, n1, n0)
  .(n_gene = length(W[[.BY$variant]]), mean_auc = mean(auc),
    pooled_auc = p[1], lo = p[2], hi = p[3], I2 = p[4]) }, by = variant]
setorder(R, -pooled_auc)
cat("\n===== 基因级方案对照 =====\n"); print(R, digits = 3)
cat("\n===== 逐队列 AUC =====\n"); print(dcast(D, variant ~ cohort, value.var = "auc"), digits = 3)
cat("\n===== 基因覆盖率 =====\n"); print(dcast(D, variant ~ cohort, value.var = "coverage"), digits = 2)
saveRDS(list(W = W, R = R, D = D, immune_genes = immune_genes), file.path(RDIR, "exp_ias2_results.rds"))
