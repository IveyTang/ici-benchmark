## ============================================================
## exp_ias4.R — 跨队列留一（LOCO）派生 ICI 应答签名
##   设计：直接在 ICI 队列上训练（免疫基因池 + LASSO/弹性网），
##         被留出的队列完全独立 —— 得到无偏的跨队列泛化估计
##   对照：固定 IFN-γ 6 基因签名
## ============================================================
suppressPackageStartupMessages({ library(data.table); library(pROC); library(metafor); library(glmnet) })
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))

coh <- readRDS(file.path(RDIR, "ici_cohorts.rds"))
zc <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }
ZC <- lapply(names(coh), function(nm) zc(coh[[nm]]$X)); names(ZC) <- names(coh)

## ---------- 候选基因：免疫基因池 ∩ 至少 5 个队列可测 ----------
read_gmt <- function(f) { ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  p <- strsplit(ln, "\t", fixed = TRUE)
  gs <- lapply(p, function(x) unique(x[-c(1, 2)])); names(gs) <- vapply(p, function(x) x[1], ""); gs }
GH <- read_gmt(file.path(DL, "msigdb/h.all.v2024.1.Hs.symbols.gmt"))
C7 <- read_gmt(file.path(DL, "msigdb/c7.all.v2024.1.Hs.symbols.gmt"))
kw <- paste(c("CD8","CD4","TREG","NK.CELL","BCELL","B.CELL","DENDRITIC","MACROPHAGE","MONOCYTE",
  "NEUTROPHIL","MAST.CELL","EOSINOPHIL","IFNG","INTERFERON","CYTOLYTIC","TCELL","T.CELL","IMMUNE",
  "GRANZYME","GZMA","GZMB","PD1","PDL1","CHECKPOINT","TH1","TH2","TH17","PLASMA.CELL","MYELOID"),
  collapse = "|")
g7 <- C7[grepl(kw, names(C7), ignore.case = TRUE)]; g7 <- g7[lengths(g7) >= 10 & lengths(g7) <= 500]
g7 <- g7[order(-lengths(g7))][1:60]
pool0 <- unique(c(unlist(GH), unlist(g7)))
npresent <- sapply(pool0, function(g) sum(sapply(ZC, function(z) g %in% rownames(z))))
pool <- names(npresent)[npresent >= 5]
pool <- setdiff(pool, c(""))                       # 去掉空名
cat("候选基因池:", length(pool), "（≥5 队列可测）\n")
cat("各队列覆盖:", paste(names(coh), sapply(ZC, function(z) sum(pool %in% rownames(z))), sep=":", collapse="  "), "\n")

## ---------- 统一到共同基因空间 ----------
P <- pool

fit_lasso <- function(train_names, alpha = 1, s = "lambda.1se", seed = 1) {
  Xl <- list(); yl <- list()
  for (nm in train_names) {
    g <- intersect(P, rownames(ZC[[nm]])); y <- coh[[nm]]$response
    keep <- !is.na(y)
    if (sum(keep) < 10) next
    Xl[[nm]] <- ZC[[nm]][g, keep, drop = FALSE]; yl[[nm]] <- y[keep]
  }
  G <- Reduce(intersect, lapply(Xl, rownames))
  X <- do.call(rbind, lapply(Xl, function(m) t(m[G, , drop = FALSE])))
  y <- unlist(yl)
  set.seed(seed)
  cv <- cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = 10,
                  standardize = FALSE)
  cc <- as.matrix(coef(cv, s = s))
  v <- cc[cc[, 1] != 0 & rownames(cc) != "(Intercept)", 1]
  list(beta = v, n = length(y), n_event = sum(y == 1), lambda = cv$lambda.min, G = G)
}
score_on <- function(beta, nm) {
  g <- intersect(names(beta), rownames(ZC[[nm]]))
  if (length(g) < 3) return(NULL)
  b <- beta[g] / sum(abs(beta[g]))
  list(score = as.numeric(t(ZC[[nm]][g, , drop = FALSE]) %*% b), n_used = length(g))
}
auc_of <- function(sc, nm) {
  y <- coh[[nm]]$response; keep <- !is.na(y) & is.finite(sc)
  if (sum(keep) < 10 || length(unique(y[keep])) < 2) return(NULL)
  r <- pROC::roc(y[keep], sc[keep], quiet = TRUE)
  data.table(cohort = nm, n = sum(keep), n_resp = sum(y[keep] == 1),
             auc = as.numeric(pROC::auc(r)),
             lo = as.numeric(pROC::ci.auc(r))[1], hi = as.numeric(pROC::ci.auc(r))[3])
}
pooled <- function(a, n1, n0) {
  d <- data.table(a = a, n1 = n1, n0 = n0); d[, logit := log(a/(1-a))]
  hm <- function(a,n1,n0){q1<-a/(2-a);q2<-2*a^2/(1+a);(a*(1-a)+(n1-1)*(q1-a^2)+(n0-1)*(q2-a^2))/(n1*n0)}
  d[, vi := hm(a, n1, n0)]; m <- rma(yi = logit, vi = vi, data = d, method = "DL")
  c(auc = plogis(m$b[1,1]), lo = plogis(m$ci.lb), hi = plogis(m$ci.ub), I2 = m$I2, p = m$pval)
}

## ---------- A. LOCO：训练于其余队列，测试被留出队列 ----------
loco_lasso <- rbindlist(lapply(names(coh), function(k) {
  f <- fit_lasso(setdiff(names(coh), k), alpha = 1, s = "lambda.1se", seed = 20261012)
  sc <- score_on(f$beta, k); if (is.null(sc)) return(NULL)
  r <- auc_of(sc$score, k); if (is.null(r)) return(NULL)
  r[, `:=`(scheme = "LOCO_LASSO_1se", n_gene = length(f$beta), n_gene_used = sc$n_used,
           train_n = f$n)]
  r
}))
loco_en <- rbindlist(lapply(names(coh), function(k) {
  f <- fit_lasso(setdiff(names(coh), k), alpha = 0.5, s = "lambda.1se", seed = 20261012)
  sc <- score_on(f$beta, k); if (is.null(sc)) return(NULL)
  r <- auc_of(sc$score, k); if (is.null(r)) return(NULL)
  r[, `:=`(scheme = "LOCO_EN_a0.5", n_gene = length(f$beta), n_gene_used = sc$n_used, train_n = f$n)]
  r
}))
## ---------- B. 固定参考签名 ----------
refs <- list(IFNG6 = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"),
             AYERS = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","HLA-DQA1","HLA-DRB1",
                       "HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT","TAP1","CXCL10"),
             CYT2 = c("GZMA","PRF1"))
loco_ref <- rbindlist(lapply(names(refs), function(rn) {
  b <- setNames(rep(1, length(intersect(refs[[rn]], P))), intersect(refs[[rn]], P))
  rbindlist(lapply(names(coh), function(k) {
    sc <- score_on(b, k); if (is.null(sc)) return(NULL)
    r <- auc_of(sc$score, k); if (is.null(r)) return(NULL)
    r[, `:=`(scheme = paste0("REF_", rn), n_gene = length(b), n_gene_used = sc$n_used, train_n = NA_integer_)]
    r
  }))
}))
## ---------- C. 仅在 IMvigor210 训练 ----------
f_imv <- fit_lasso("IMvigor210", alpha = 1, s = "lambda.1se", seed = 20261012)
imv_only <- rbindlist(lapply(setdiff(names(coh), "IMvigor210"), function(k) {
  sc <- score_on(f_imv$beta, k); if (is.null(sc)) return(NULL)
  r <- auc_of(sc$score, k); if (is.null(r)) return(NULL)
  r[, `:=`(scheme = "TrainIMvigor210", n_gene = length(f_imv$beta), n_gene_used = sc$n_used,
           train_n = f_imv$n)]
  r
}))

ALL <- rbindlist(list(loco_lasso, loco_en, loco_ref, imv_only), fill = TRUE)
S <- ALL[, { p <- pooled(auc, n_resp, n - n_resp)
  .(n_gene = n_gene[1], k = .N, mean_auc = mean(auc), pooled_auc = p[1],
    lo = p[2], hi = p[3], I2 = p[4], pval = p[5]) }, by = scheme]
setorder(S, -pooled_auc)
cat("\n===== 方案对比（留出队列）=====\n"); print(S, digits = 3)
cat("\n===== 逐队列 AUC =====\n"); print(dcast(ALL, scheme ~ cohort, value.var = "auc"), digits = 3)
cat("\n===== IMvigor210 训练模型选中的基因（", length(f_imv$beta), " 个）=====\n")
cat(paste(names(f_imv$beta), collapse = ", "), "\n")
saveRDS(list(ALL = ALL, S = S, f_imv = f_imv, pool = pool), file.path(RDIR, "exp_ias4_results.rds"))
