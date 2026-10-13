## ============================================================
## exp_ias3.R — 以 IFN-γ / T 细胞炎症轴为锚点，在 TCGA 中派生紧凑 IAS
##   动机：GSVA 表型分离弱(0.09 SD)，监督权重不泛化；改用生物学锚点做标签
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
ZC <- lapply(names(coh), function(nm) zc(coh[[nm]]$X)); names(ZC) <- names(coh)

## ---------- 签名库 ----------
read_gmt <- function(f) { ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  p <- strsplit(ln, "\t", fixed = TRUE)
  gs <- lapply(p, function(x) unique(x[-c(1, 2)])); names(gs) <- vapply(p, function(x) x[1], ""); gs }
GH <- read_gmt(file.path(DL, "msigdb/h.all.v2024.1.Hs.symbols.gmt"))
IFNG_HALL <- GH[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]
IFNG_6 <- c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")
AYERS <- c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","HLA-DQA1","HLA-DRB1","HLA-E",
           "IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT","TAP1","CXCL10")

## ---------- TCGA 锚点分数 ----------
anchor_genes <- intersect(unique(c(IFNG_HALL, IFNG_6, AYERS)), rownames(TX))
ZA <- zc(TX[anchor_genes, ])
anchor <- colMeans(ZA[intersect(IFNG_HALL, rownames(ZA)), ])   # 用 Hallmark IFN-γ 全集作锚
names(anchor) <- colnames(TX)
cat("锚点基因:", length(intersect(IFNG_HALL, rownames(TX))), " IFN-γ score sd =", round(sd(anchor),3), "\n")
## TCGA 免疫基因池
Gpool <- intersect(unique(c(IFNG_HALL, AYERS,
  unlist(GH[c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_ALLOGRAFT_REJECTION",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_COMPLEMENT")]))), rownames(TX))
cat("候选基因池:", length(Gpool), "\n")

yq <- as.integer(anchor > quantile(anchor, .75)); yq[anchor < quantile(anchor, .25)] <- 0L
names(yq) <- names(anchor); ok <- !is.na(yq) & (anchor > quantile(anchor,.75) | anchor < quantile(anchor,.25))
Xf <- t(zc(TX[Gpool, ]))

W <- list()
for (al in c(0.5, 1)) for (s in c("lambda.1se","lambda.min")) {
  set.seed(20261012)
  cv <- cv.glmnet(Xf[ok, ], yq[ok], family = "binomial", alpha = al, nfolds = 10)
  cc <- as.matrix(coef(cv, s = s)); v <- cc[cc[,1] != 0 & rownames(cc) != "(Intercept)", 1]
  W[[sprintf("EN_a%.1f_%s", al, sub("lambda\\.","",s))]] <- v
}
## 相关性加权 top-N
cr <- as.numeric(cor(Xf, anchor)); names(cr) <- colnames(Xf)
for (N in c(10, 15, 20, 30)) W[[paste0("COR_top", N)]] <- cr[order(-abs(cr))][1:N]
## 参考签名
W[["REF_IFNG6"]] <- setNames(rep(1, length(intersect(IFNG_6, rownames(TX)))),
                             intersect(IFNG_6, rownames(TX)))
W[["REF_AYERS"]] <- setNames(rep(1, length(intersect(AYERS, rownames(TX)))),
                             intersect(AYERS, rownames(TX)))
W[["REF_IFNG_HALL"]] <- setNames(rep(1, length(intersect(IFNG_HALL, rownames(TX)))),
                                 intersect(IFNG_HALL, rownames(TX)))

pooled <- function(a, n1, n0) {
  d <- data.table(a = a, n1 = n1, n0 = n0); d[, logit := log(a/(1-a))]
  hm <- function(a,n1,n0){q1<-a/(2-a);q2<-2*a^2/(1+a);(a*(1-a)+(n1-1)*(q1-a^2)+(n0-1)*(q2-a^2))/(n1*n0)}
  d[, vi := hm(a, n1, n0)]; m <- rma(yi = logit, vi = vi, data = d, method = "DL")
  c(auc = plogis(m$b[1,1]), lo = plogis(m$ci.lb), hi = plogis(m$ci.ub), I2 = m$I2, p = m$pval)
}
detail <- function(bv, label) rbindlist(lapply(names(coh), function(nm) {
  g <- intersect(names(bv), rownames(ZC[[nm]])); if (length(g) < 3) return(NULL)
  b <- bv[g] / sum(abs(bv[g]))
  s <- as.numeric(t(ZC[[nm]][g, , drop = FALSE]) %*% b)
  yv <- coh[[nm]]$response; keep <- !is.na(yv) & is.finite(s)
  if (sum(keep) < 10 || length(unique(yv[keep])) < 2) return(NULL)
  data.table(variant = label, cohort = nm, auc = as.numeric(pROC::auc(pROC::roc(yv[keep], s[keep], quiet = TRUE))),
             n1 = sum(yv[keep]==1), n0 = sum(yv[keep]==0))
}))
D <- rbindlist(lapply(names(W), function(nm) detail(W[[nm]], nm)))
R <- D[, { p <- pooled(auc, n1, n0)
  .(n_gene = length(W[[.BY$variant]]), mean_auc = mean(auc), pooled_auc = p[1],
    lo = p[2], hi = p[3], I2 = p[4], pval = p[5]) }, by = variant]
setorder(R, -pooled_auc)
cat("\n===== 结果 =====\n"); print(R, digits = 3)
cat("\n===== 逐队列 =====\n"); print(dcast(D, variant ~ cohort, value.var = "auc"), digits = 3)
saveRDS(list(W = W, R = R, D = D, anchor_genes = anchor_genes), file.path(RDIR, "exp_ias3_results.rds"))
