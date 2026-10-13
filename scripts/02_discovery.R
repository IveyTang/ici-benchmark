## ============================================================
## 02_discovery.R — TCGA-HNSC 发现：免疫表型、生存关联、IAS 构建
##   * 富集打分统一使用 GSVA（核密度法）——ssGSEA 会压缩样本间结构
##   * 免疫表型：对 GSVA 分数矩阵做一致性聚类（Spearman 距离, k = 2）
##   * IAS：LASSO 逻辑回归，仅在 70% 训练集上拟合，报告 30% 测试集 AUC（避免循环论证）
## ============================================================
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(msigdbr); library(GSVA); library(ConsensusClusterPlus)
  library(survival); library(glmnet); library(pROC)
})

tc <- readRDS(file.path(RDIR, "tcga_hnsc.rds"))
X <- tc$X; surv <- tc$surv; clin <- tc$clin; tmb <- tc$tmb
msg("TCGA 输入: ", nrow(X), " genes x ", ncol(X), " samples")

## ---------- 1. 构建免疫签名库（本地 MSigDB GMT，避免 zenodo 依赖）----------
msg("构建免疫签名库（MSigDB 2024.1.Hs 本地 GMT）")
read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE)
  ln <- ln[nzchar(ln)]
  parts <- strsplit(ln, "\t", fixed = TRUE)
  gs <- lapply(parts, function(x) unique(x[-c(1, 2)]))
  names(gs) <- vapply(parts, function(x) x[1], "")
  gs[lengths(gs) > 0]
}
gsc_h <- read_gmt(file.path(DL, "msigdb/h.all.v2024.1.Hs.symbols.gmt"))
c7_all <- read_gmt(file.path(DL, "msigdb/c7.all.v2024.1.Hs.symbols.gmt"))
kw <- paste(c("CD8", "CD4", "TREG", "REGULATORY.T", "NK.CELL", "NKT", "BCELL", "B.CELL",
              "DENDRITIC", "MACROPHAGE", "MONOCYTE", "NEUTROPHIL", "MAST.CELL", "EOSINOPHIL",
              "IFNG", "INTERFERON", "CYTOLYTIC", "TCELL", "T.CELL", "IMMUNE",
              "GRANZYME", "GZMA", "GZMB", "PD1", "PDL1", "CHECKPOINT", "TH1", "TH2", "TH17",
              "CENTRAL.MEMORY", "EFFECTOR.MEMORY", "PLASMA.CELL", "MYELOID"), collapse = "|")
gsc_c7 <- c7_all[grepl(kw, names(c7_all), ignore.case = TRUE)]
gsc_c7 <- gsc_c7[lengths(gsc_c7) >= 10 & lengths(gsc_c7) <= 500]
if (length(gsc_c7) > 60) gsc_c7 <- gsc_c7[order(-lengths(gsc_c7))][1:60]
msg("  Hallmark 集合 ", length(gsc_h), " 个; C7 免疫相关 ", length(gsc_c7), " 个")
gsc_all <- c(gsc_h, gsc_c7)
gsc_all <- lapply(gsc_all, unique)
gsc_all <- gsc_all[lengths(gsc_all) >= 10]

## 限制在 TCGA 检出基因内，保证打分稳定
gsc <- lapply(gsc_all, function(g) intersect(g, rownames(X)))
gsc <- gsc[lengths(gsc) >= 10]
msg("  最终参与评分的签名: ", length(gsc), " 个（中位基因数 ", median(lengths(gsc)), "）")

## ---------- 2. GSVA 富集打分 ----------
f_gs <- file.path(RDIR, "tcga_enrich_gsva.rds")
if (file.exists(f_gs)) {
  msg("加载缓存的 GSVA 分数矩阵")
  gs <- readRDS(f_gs)
} else {
  msg("计算 GSVA 富集分数（核密度法）")
  gs <- run_enrich(X, gsc, method = "gsva")
  saveRDS(gs, file.path(RDIR, "tcga_enrich_gsva.rds"))
}
msg("  分数矩阵: ", nrow(gs), " 签名 x ", ncol(gs), " 样本")
zg <- t(scale(t(gs))); zg[!is.finite(zg)] <- 0      # 每个签名在样本间 z 化
imm_score <- colMeans(zg)                            # 全局免疫活化分数

## ---------- 3. 一致性聚类（k = 2）----------
msg("一致性聚类（k = 2-6, Spearman 距离, 1000 次重采样）")
set.seed(20261012)
cc <- ConsensusClusterPlus(as.matrix(gs), maxK = 6, reps = 1000, pItem = 0.8,
                           pFeature = 1, clusterAlg = "hc", distance = "spearman",
                           innerLinkage = "average", finalLinkage = "average",
                           seed = 20261012, plot = NULL, verbose = FALSE)
## CDF 与 delta-area（聚类数目选择的定量依据）
cdf_at <- function(k, x) mean(cc[[k]]$consensusMatrix[lower.tri(cc[[k]]$consensusMatrix)] <= x)
grid <- seq(0.01, 0.99, by = 0.01)
cdf_k <- sapply(2:6, function(k) sapply(grid, function(x) cdf_at(k, x)))
areas <- apply(cdf_k, 2, function(v) sum(diff(grid) * (head(v, -1) + tail(v, -1)) / 2))
delta <- c(NA, -diff(areas))
png(file.path(FDIR, "fig0_consensus_cdf.png"), 2800, 1200, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
matplot(grid, cdf_k, type = "l", lty = 1, lwd = 2, col = 1:5,
        xlab = "Consensus index", ylab = "CDF", main = "Consensus CDF")
legend("bottomright", paste0("k = ", 2:6), col = 1:5, lwd = 2, bty = "n", cex = 0.8)
plot(2:6, delta, type = "b", pch = 19, lwd = 2, xlab = "k", ylab = "Delta area",
     main = "Relative change in area under CDF")
dev.off()
msg("  delta area: ", paste(2:6, signif(delta, 4), sep = "=", collapse = " "))
msg("  各 k 簇规模: ", paste(sapply(2:6, function(k)
  paste0("k", k, "[", paste(as.integer(table(cc[[k]]$consensusClass)), collapse = "/"), "]")),
  collapse = " "))

clu <- cc[[2]]$consensusClass
names(clu) <- colnames(gs)
mu <- tapply(imm_score, clu, mean); hi <- names(mu)[which.max(mu)]
pheno <- as.integer(clu == hi); names(pheno) <- names(clu)      # 1 = immune-activated
msg("  k=2 簇规模: ", paste(names(table(clu)), as.integer(table(clu)), sep = "=", collapse = " "),
    "；两簇免疫分数均值 ", paste(round(mu, 3), collapse = " vs "),
    "；immune-activated = ", sum(pheno), " (", round(100 * mean(pheno), 1), "%)")
saveRDS(list(clu = clu, pheno = pheno, imm_score = imm_score, delta = delta),
        file.path(RDIR, "tcga_phenotype.rds"))

## ---------- 4. 生存关联 ----------
msg("生存分析")
sd <- merge(data.table(sample = names(pheno), pheno = pheno),
            as.data.table(surv), by = "sample")
sd <- sd[!is.na(OS.time) & OS.time > 0]
km <- survdiff(Surv(OS.time, OS) ~ pheno, data = sd)
fit <- coxph(Surv(OS.time, OS) ~ pheno, data = sd)
s_cox <- summary(fit)
res_surv <- data.table(
  test = c("log-rank", "Cox HR (immune-activated vs cold)"),
  stat = c(km$chisq, s_cox$coefficients[1, "exp(coef)"]),
  lo = c(NA, s_cox$conf.int[1, "lower .95"]), hi = c(NA, s_cox$conf.int[1, "upper .95"]),
  p = c(1 - pchisq(km$chisq, 1), s_cox$coefficients[1, "Pr(>|z|)"]))
msg("  log-rank P = ", signif(res_surv$p[1], 3),
    "; HR = ", signif(res_surv$stat[2], 3),
    " (", signif(res_surv$lo[2], 3), "-", signif(res_surv$hi[2], 3), "), P = ",
    signif(res_surv$p[2], 3))

## 多因素校正（年龄、性别、分期）
cl <- as.data.table(clin)
cl[, age := as.numeric(age_at_index.demographic)]
cl[, sex := gender.demographic]
cl[, stage_iv := grepl("IV", ajcc_pathologic_stage.diagnoses)]
mv <- merge(sd, cl[, .(sample, age, sex, stage_iv)], by = "sample", all.x = TRUE)
mv <- mv[!is.na(age)]
fit2 <- coxph(Surv(OS.time, OS) ~ pheno + age + sex + stage_iv, data = mv)
s2 <- summary(fit2)
msg("  多因素校正 HR = ", signif(s2$coefficients["pheno", "exp(coef)"], 3),
    " (", signif(s2$conf.int["pheno", "lower .95"], 3), "-",
    signif(s2$conf.int["pheno", "upper .95"], 3), "), P = ",
    signif(s2$coefficients["pheno", "Pr(>|z|)"], 3), "  [n=", nrow(mv), "]")

## ---------- 5. TMB 与 PD-L1 ----------
tb <- merge(data.table(sample = names(pheno), pheno = pheno),
            as.data.table(tmb)[, .(sample, TMB)], by = "sample")
wt <- wilcox.test(TMB ~ pheno, data = tb)
msg("  TMB: immune-activated 中位 ", round(median(tb[pheno == 1]$TMB), 2),
    " vs cold ", round(median(tb[pheno == 0]$TMB), 2), " mut/Mb, P = ", signif(wt$p.value, 3))
pdl1 <- if ("CD274" %in% rownames(X)) {
  d <- data.table(sample = colnames(X), CD274 = as.numeric(X["CD274", ]))
  d <- merge(d, data.table(sample = names(pheno), pheno = pheno), by = "sample")
  w2 <- wilcox.test(CD274 ~ pheno, data = d)
  list(p = w2$p.value, hi = median(d[pheno == 1]$CD274), lo = median(d[pheno == 0]$CD274))
} else NULL
if (!is.null(pdl1)) msg("  CD274(PD-L1): ", round(pdl1$hi, 2), " vs ", round(pdl1$lo, 2),
                        ", P = ", signif(pdl1$p, 3))

## ---------- 6. LASSO 构建 IAS（训练集拟合 / 测试集评估）----------
msg("LASSO 特征选择构建 IAS（70/30 分层拆分）")
Xf <- t(zg)                                  # samples x signatures
y <- pheno[rownames(Xf)]
msg("  类别分布: ", paste(names(table(y)), as.integer(table(y)), sep = "=", collapse = " "))
stopifnot(min(table(y)) >= 15)
set.seed(20261012)
idx1 <- which(y == 1); idx0 <- which(y == 0)
tr <- c(sample(idx1, floor(0.7 * length(idx1))), sample(idx0, floor(0.7 * length(idx0))))
te <- setdiff(seq_along(y), tr)
msg("  训练 n = ", length(tr), "（应答/活化 ", sum(y[tr] == 1), "）；测试 n = ", length(te),
    "（活化 ", sum(y[te] == 1), "）")
cvf <- cv.glmnet(Xf[tr, ], y[tr], family = "binomial", alpha = 1, nfolds = 10)
coefs <- as.matrix(coef(cvf, s = "lambda.min"))
sel <- rownames(coefs)[coefs[, 1] != 0 & rownames(coefs) != "(Intercept)"]
beta <- coefs[sel, 1]
msg("  选中特征 ", length(sel), " 个 @ lambda.min = ", signif(cvf$lambda.min, 4))

## 测试集评估（非循环）
ias_te <- as.numeric(Xf[te, sel, drop = FALSE] %*% beta)
roc_te <- pROC::roc(y[te], ias_te, quiet = TRUE)
msg("  测试集 AUC = ", round(as.numeric(pROC::auc(roc_te)), 3),
    " (95% CI ", round(ci.auc(roc_te)[1], 3), "-", round(ci.auc(roc_te)[3], 3), ")")
## 重复 5 折交叉验证 AUC（训练集内部）
set.seed(7)
cvauc <- replicate(10, {
  f <- sample(rep(1:5, length.out = length(tr)))
  mean(sapply(1:5, function(k) {
    i <- which(f != k)
    m <- glmnet(Xf[tr[i], sel, drop = FALSE], y[tr[i]], family = "binomial", alpha = 1,
                lambda = cvf$lambda.min)
    s <- as.numeric(predict(m, Xf[tr[-i], sel, drop = FALSE]))
    as.numeric(pROC::auc(pROC::roc(y[tr[-i]], s, quiet = TRUE)))
  }))
})
msg("  训练集重复 CV AUC = ", round(mean(cvauc), 3), " ± ", round(sd(cvauc), 3))

## 全样本 IAS（供后续图表与生存分析使用）
ias <- as.numeric(scale(as.numeric(Xf[, sel, drop = FALSE] %*% beta)))
names(ias) <- rownames(Xf)

## IAS 与生存
sd2 <- merge(data.table(sample = names(ias), IAS = ias), as.data.table(surv), by = "sample")
fit3 <- coxph(Surv(OS.time, OS) ~ IAS, data = sd2)
s3 <- summary(fit3)
msg("  IAS 每 SD 的 OS HR = ", signif(s3$coefficients[1, "exp(coef)"], 3),
    " (", signif(s3$conf.int[1, "lower .95"], 3), "-", signif(s3$conf.int[1, "upper .95"], 3),
    "), P = ", signif(s3$coefficients[1, "Pr(>|z|)"], 3))

## ---------- 7. 保存 ----------
res <- list(
  signatures = names(beta), beta = beta,
  gene_sets = gsc[names(beta)],
  enrich_method = "gsva",
  surv = res_surv,
  surv_multi = list(hr = s2$coefficients["pheno", "exp(coef)"],
                    lo = s2$conf.int["pheno", "lower .95"], hi = s2$conf.int["pheno", "upper .95"],
                    p = s2$coefficients["pheno", "Pr(>|z|)"], n = nrow(mv)),
  tmb = list(p = wt$p.value, hi = median(tb[pheno == 1]$TMB), lo = median(tb[pheno == 0]$TMB)),
  pdl1 = pdl1,
  auc_tcga_test = as.numeric(pROC::auc(roc_te)),
  auc_tcga_test_ci = as.numeric(ci.auc(roc_te))[c(1, 3)],
  auc_tcga_cv = c(mean(cvauc), sd(cvauc)),
  n_train = length(tr), n_test = length(te),
  ias_os = list(hr = s3$coefficients[1, "exp(coef)"],
                lo = s3$conf.int[1, "lower .95"], hi = s3$conf.int[1, "upper .95"],
                p = s3$coefficients[1, "Pr(>|z|)"]),
  n_samples = ncol(X), n_genes = nrow(X), n_events = sum(surv$OS),
  n_sets_total = length(gsc), delta_area = delta,
  cluster_sizes = as.integer(table(clu)),
  followup_median_months = median(surv$OS.time) / 30.44)
saveRDS(res, file.path(RDIR, "discovery_results.rds"))

## 输出表
fwrite(data.table(sample = names(ias), IAS = ias, phenotype = y,
                  split = ifelse(seq_along(y) %in% tr, "train", "test")),
       file.path(TDIR, "tcga_ias_scores.csv"))
fwrite(data.table(signature = names(beta), beta = as.numeric(beta),
                  n_genes = lengths(gsc[names(beta)])),
       file.path(TDIR, "ias_coefficients.csv"))
fwrite(res_surv, file.path(TDIR, "tcga_survival.csv"))
msg("DISCOVERY DONE")
