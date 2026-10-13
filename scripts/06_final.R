## ============================================================
## 06_final.R — 最终权威结果集（数值 + 图 + 表）
##   * IAS（主）：TCGA-HNSC 中以致 IFN-γ 轴为锚、免疫基因池 LASSO(λ.1se) 派生
##   * 参照签名：IFN-γ 6 基因 / Ayers 18 基因 / 细胞溶解 2 基因
##   * 跨 6 个 ICI 队列：逐队列 AUC / OR / HR + 随机效应合并 + LOCO + 置换
## ============================================================
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(pROC); library(survival); library(metafor); library(glmnet); library(GSVA)
})

tc   <- readRDS(file.path(RDIR, "tcga_hnsc.rds"))
disc <- readRDS(file.path(RDIR, "discovery_results.rds"))
phe  <- readRDS(file.path(RDIR, "tcga_phenotype.rds"))
coh  <- readRDS(file.path(RDIR, "ici_cohorts.rds"))
zc <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }

## ---------- 0. 引用签名 ----------
REF <- list(
  IFNG6 = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"),
  AYERS = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","HLA-DQA1","HLA-DRB1","HLA-E",
            "IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT","TAP1","CXCL10"),
  CYT2  = c("GZMA","PRF1"))

## ---------- 1. TCGA 中派生 IAS（IFN-γ 轴锚点 + LASSO）----------
read_gmt <- function(f) { ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  p <- strsplit(ln, "\t", fixed = TRUE)
  gs <- lapply(p, function(x) unique(x[-c(1, 2)])); names(gs) <- vapply(p, function(x) x[1], ""); gs }
GH <- read_gmt(file.path(DL, "msigdb/h.all.v2024.1.Hs.symbols.gmt"))
IFNG_HALL <- intersect(GH[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], rownames(tc$X))
Gpool <- intersect(unique(c(IFNG_HALL, REF$AYERS,
  unlist(GH[c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_ALLOGRAFT_REJECTION",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_COMPLEMENT")]))), rownames(tc$X))
ZP <- zc(tc$X[Gpool, ])
anchor <- colMeans(zc(tc$X[IFNG_HALL, ])); names(anchor) <- colnames(tc$X)
yq <- as.integer(anchor > quantile(anchor, .75)); names(yq) <- names(anchor)
ok <- anchor > quantile(anchor, .75) | anchor < quantile(anchor, .25)
set.seed(20261012)
cv <- cv.glmnet(t(ZP)[ok, ], yq[ok], family = "binomial", alpha = 1, nfolds = 10)
cc <- as.matrix(coef(cv, s = "lambda.1se"))
beta <- cc[cc[, 1] != 0 & rownames(cc) != "(Intercept)", 1]
msg("IAS = ", length(beta), " 个基因（IFN-γ 轴锚点, LASSO λ.1se = ", signif(cv$lambda.1se, 4), "）")
Zt <- zc(tc$X[names(beta), ])
ias_tc <- as.numeric(scale(as.numeric(t(Zt) %*% (beta / sum(abs(beta))))))
names(ias_tc) <- colnames(tc$X)

## ---------- 2. 各队列打分 ----------
scorer <- function(bv, nm) {
  Z <- zc(coh[[nm]]$X); g <- intersect(names(bv), rownames(Z))
  if (length(g) < 3) return(NULL)
  b <- bv[g] / sum(abs(bv[g]))
  s <- as.numeric(t(Z[g, , drop = FALSE]) %*% b)
  list(score = as.numeric(scale(s)), n_used = length(g), coverage = length(g) / length(bv))
}
CAND <- c(list(IAS = beta), lapply(REF, function(g) setNames(rep(1, length(g)), g)))

rows <- list()
for (cn in names(CAND)) for (nm in names(coh)) {
  r <- scorer(CAND[[cn]], nm); if (is.null(r)) next
  o <- coh[[nm]]; y <- o$response
  keep <- !is.na(y) & is.finite(r$score)
  if (sum(keep) < 10 || length(unique(y[keep])) < 2) next
  rr <- pROC::roc(y[keep], r$score[keep], quiet = TRUE)
  a <- as.numeric(pROC::auc(rr)); ci <- as.numeric(pROC::ci.auc(rr))
  gl <- glm(y[keep] ~ r$score[keep], family = binomial)
  or <- exp(coef(gl)[2]); orl <- exp(confint(gl)[2, 1]); orh <- exp(confint(gl)[2, 2])
  hr <- hrl <- hrh <- hrp <- NA_real_
  if (!is.null(o$os_time) && !is.null(o$os_event)) {
    d <- data.table(t = o$os_time[keep], e = o$os_event[keep], s = r$score[keep])
    d <- d[!is.na(t) & !is.na(e) & t > 0]
    if (nrow(d) > 20 && sum(d$e) > 5) {
      f <- tryCatch(coxph(Surv(t, e) ~ s, data = d), error = function(e) NULL)
      if (!is.null(f)) { sf <- summary(f); hr <- sf$coefficients[1, "exp(coef)"]
        hrl <- sf$conf.int[1, "lower .95"]; hrh <- sf$conf.int[1, "upper .95"]
        hrp <- sf$coefficients[1, "Pr(>|z|)"] }
    }
  }
  rows[[length(rows) + 1]] <- data.table(
    signature = cn, cohort = nm, cancer = o$cancer, platform = o$platform,
    n = sum(keep), n_resp = sum(y[keep] == 1), n_nonresp = sum(y[keep] == 0),
    orr = mean(y[keep] == 1), n_gene_used = r$n_used, coverage = r$coverage,
    auc = a, auc_lo = ci[1], auc_hi = ci[3], or = or, or_lo = orl, or_hi = orh,
    os_hr = hr, os_hr_lo = hrl, os_hr_hi = hrh, os_p = hrp)
}
perf <- rbindlist(rows)
fwrite(perf, file.path(TDIR, "final_cohort_performance.csv"))
cat("\n===== 逐队列表现 =====\n"); print(perf[, .(signature, cohort, n, n_resp, auc = round(auc,3),
  ci = sprintf("%.3f-%.3f", auc_lo, auc_hi), or = round(or,2))])

## ---------- 3. 随机效应合并 ----------
hanley_var <- function(auc, n1, n0) { q1 <- auc/(2-auc); q2 <- 2*auc^2/(1+auc)
  (auc*(1-auc) + (n1-1)*(q1-auc^2) + (n0-1)*(q2-auc^2))/(n1*n0) }
pool_all <- rbindlist(lapply(unique(perf$signature), function(cn) {
  d <- copy(perf[signature == cn]); d[, logit := log(auc/(1-auc))]
  d[, vi := hanley_var(auc, n_resp, n_nonresp)]
  m <- rma(yi = logit, vi = vi, data = d, method = "DL")
  d[, vi_or := (log(or_hi) - log(or_lo))^2 / (2*1.96)^2]
  mo <- rma(yi = log(or), vi = vi_or, data = d, method = "DL")
  data.table(signature = cn, k = nrow(d), total_n = sum(d$n),
    pooled_auc = plogis(m$b[1,1]), lo = plogis(m$ci.lb), hi = plogis(m$ci.ub),
    I2 = m$I2, Q = m$QE, Qp = m$QEp, tau2 = m$tau2,
    pooled_or = exp(mo$b[1,1]), or_lo = exp(mo$ci.lb), or_hi = exp(mo$ci.ub), or_p = mo$pval)
}))
setorder(pool_all, -pooled_auc)
fwrite(pool_all, file.path(TDIR, "final_pooled.csv"))
cat("\n===== 合并效应量 =====\n"); print(pool_all, digits = 3)

## ---------- 4. LOCO（留一队列）----------
loco <- rbindlist(lapply(unique(perf$signature), function(cn) {
  d <- copy(perf[signature == cn]); d[, logit := log(auc/(1-auc))]
  d[, vi := hanley_var(auc, n_resp, n_nonresp)]
  rbindlist(lapply(unique(d$cohort), function(k) {
    m <- tryCatch(rma(yi = logit, vi = vi, data = d[cohort != k], method = "DL"), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.table(signature = cn, excluded = k, pooled_auc = plogis(m$b[1,1]),
               lo = plogis(m$ci.lb), hi = plogis(m$ci.ub), I2 = m$I2)
  }))
}))
fwrite(loco, file.path(TDIR, "final_leave_one_out.csv"))
cat("\n===== LOCO 范围 =====\n")
print(loco[, .(min = round(min(pooled_auc),3), max = round(max(pooled_auc),3)), by = signature])

## ---------- 5. 置换检验（IAS）----------
msg("置换检验（IAS, 200 次）")
tcX <- tc$X[intersect(rownames(tc$X), names(beta)), , drop = FALSE]
tz <- tcX
y <- phe$pheno[colnames(tz)]
obs <- as.numeric(scale(as.numeric(t(zc(tz)) %*% (beta / sum(abs(beta))))))
auc_obs <- as.numeric(pROC::auc(pROC::roc(y, obs, quiet = TRUE)))
set.seed(20261012); nperm <- 200; pa <- rep(NA_real_, nperm)
for (i in seq_len(nperm)) {
  bg <- setNames(rnorm(length(beta)), sample(rownames(tz), length(beta), replace = FALSE))
  s <- as.numeric(scale(as.numeric(t(zc(tz)) %*% (bg / sum(abs(bg))))))
  pa[i] <- as.numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE)))
}
p_perm <- mean(pa >= auc_obs)
msg("  IAS 对 TCGA 免疫表型的 AUC = ", round(auc_obs, 3), "；置换中位 = ", round(median(pa), 3),
    "；经验 P = ", signif(max(p_perm, 1/nperm), 3))

## ---------- 6. 图 ----------
pcol <- colorRampPalette(c("#2c5aa0", "white", "#c0392b"))(100)
## 图1 免疫表型热图（GSVA）
gs <- readRDS(file.path(RDIR, "tcga_enrich_gsva.rds"))
o <- order(phe$pheno, phe$imm_score)
png(file.path(FDIR, "fig1_phenotype_heatmap.png"), 3400, 2300, res = 300)
par(mar = c(4, 13, 3, 5))
zl <- t(scale(t(gs[, o]))); zl[zl > 2] <- 2; zl[zl < -2] <- -2
image(1:ncol(zl), 1:nrow(zl), t(zl[nrow(zl):1, ]), col = pcol, axes = FALSE,
      xlab = "TCGA-HNSC primary tumours (n = 520)", ylab = "",
      main = "Immune gene-set enrichment scores (GSVA, z-scored)")
axis(2, at = 1:nrow(zl), labels = rev(rownames(zl)), las = 1, cex.axis = 0.4, tick = FALSE)
abline(v = sum(phe$pheno[o] == 0) + 0.5, lwd = 2.5)
legend("topright", c("Immune-cold", "Immune-activated"), lwd = 2.5, bty = "n", cex = 0.85)
dev.off()
## 图2 KM
sd <- merge(data.table(sample = names(phe$pheno), pheno = phe$pheno), as.data.table(tc$surv), by = "sample")
sd <- sd[!is.na(OS.time) & OS.time > 0]; sd[, m := OS.time / 30.44]
png(file.path(FDIR, "fig2_km_phenotype.png"), 2400, 2100, res = 300)
fit <- survfit(Surv(m, OS) ~ pheno, data = sd)
plot(fit, col = c("#2c5aa0", "#c0392b"), lwd = 2.5, xlab = "Time (months)",
     ylab = "Overall survival", mark.time = TRUE,
     main = sprintf("TCGA-HNSC overall survival by immune phenotype (log-rank P = %.3g)",
                    1 - pchisq(survdiff(Surv(m, OS) ~ pheno, data = sd)$chisq, 1)))
legend("topright", c(sprintf("Immune-cold (n = %d)", sum(sd$pheno == 0)),
                     sprintf("Immune-activated (n = %d)", sum(sd$pheno == 1))),
       col = c("#2c5aa0", "#c0392b"), lwd = 2.5, bty = "n")
dev.off()
## 图3 TMB / PD-L1
tb <- merge(data.table(sample = names(phe$pheno), pheno = phe$pheno),
            as.data.table(tc$tmb)[, .(sample, TMB)], by = "sample")
png(file.path(FDIR, "fig3_tmb_pdl1.png"), 2600, 1300, res = 300)
par(mfrow = c(1, 2))
boxplot(TMB ~ pheno, data = tb, col = c("#cfd8e3", "#f2c4c0"), names = c("Cold", "Activated"),
        ylab = "TMB (mut/Mb)", main = "Tumour mutational burden",
        xlab = sprintf("Wilcoxon P = %.3g", wilcox.test(TMB ~ pheno, data = tb)$p.value))
d <- data.table(CD274 = as.numeric(tc$X["CD274", ]), pheno = phe$pheno)
boxplot(CD274 ~ pheno, data = d, col = c("#cfd8e3", "#f2c4c0"), names = c("Cold", "Activated"),
        ylab = "CD274 (PD-L1), log2(FPKM+1)", main = "PD-L1 expression",
        xlab = sprintf("Wilcoxon P = %.3g", wilcox.test(CD274 ~ pheno, data = d)$p.value))
dev.off()
## 图4 ROC（IAS + IFN-γ 6）
sc <- rbindlist(lapply(names(coh), function(nm) {
  r1 <- scorer(CAND$IAS, nm); r2 <- scorer(CAND$IFNG6, nm)
  if (is.null(r1)) return(NULL)
  data.table(cohort = nm, id = colnames(coh[[nm]]$X), IAS = r1$score,
             IFNG6 = if (is.null(r2)) NA_real_ else r2$score,
             response = coh[[nm]]$response)
}))
fwrite(sc, file.path(TDIR, "final_ias_scores.csv"))
png(file.path(FDIR, "fig4_roc_by_cohort.png"), 2700, 2100, res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (nm in names(coh)) {
  dd <- sc[cohort == nm]; k <- !is.na(dd$response)
  r1 <- pROC::roc(dd$response[k], dd$IAS[k], quiet = TRUE)
  r2 <- pROC::roc(dd$response[k], dd$IFNG6[k], quiet = TRUE)
  plot(r1, col = "#c0392b", lwd = 2, legacy.axes = TRUE,
       main = sprintf("%s (n=%d)", nm, sum(k)))
  lines(r2, col = "#2c5aa0", lwd = 2)
  legend("bottomright", c(sprintf("IAS (AUC=%.2f)", as.numeric(pROC::auc(r1))),
                          sprintf("IFN-\u03b3 (AUC=%.2f)", as.numeric(pROC::auc(r2)))),
         col = c("#c0392b", "#2c5aa0"), lwd = 2, bty = "n", cex = 0.8)
}
dev.off()
## 图5 森林图
png(file.path(FDIR, "fig5_forest_auc.png"), 2700, 1800, res = 300)
par(mar = c(4.5, 13, 3, 4))
dd <- perf[signature == "IAS"][order(auc)]
plot(dd$auc, nrow(dd):1, xlim = c(0.3, 0.95), pch = 18, cex = 1.9, xlab = "AUC (95% CI)",
     ylab = "", yaxt = "n", main = "IAS discrimination of ICI response by cohort")
segments(dd$auc_lo, nrow(dd):1, dd$auc_hi, nrow(dd):1, lwd = 2)
abline(v = 0.5, lty = 3, col = "grey50")
pa1 <- pool_all[signature == "IAS"]
abline(v = pa1$pooled_auc, lty = 2, col = "#c0392b", lwd = 2)
axis(2, at = nrow(dd):1, labels = sprintf("%s (%s, n=%d)", dd$cohort, dd$cancer, dd$n),
     las = 1, cex.axis = 0.85)
mtext(sprintf("Pooled AUC = %.3f (%.3f-%.3f), I\u00b2 = %.0f%%", pa1$pooled_auc, pa1$lo, pa1$hi, pa1$I2),
      side = 1, line = 3, cex = 0.9)
dev.off()

## ---------- 7. 表 ----------
cl <- as.data.table(tc$clin)
cl[, age := as.numeric(age_at_index.demographic)]
cl[, stage_iv := grepl("IV", ajcc_pathologic_stage.diagnoses)]
cl <- merge(cl, data.table(sample = names(phe$pheno), pheno = phe$pheno), by = "sample")
cl <- merge(cl, as.data.table(tc$tmb)[, .(sample, TMB)], by = "sample", all.x = TRUE)
cl <- merge(cl, as.data.table(tc$surv), by = "sample", all.x = TRUE)
g <- function(sub, col) sub[[col]]
mk <- function(sub) c(
  n = as.character(nrow(sub)),
  age = sprintf("%.0f (%.0f-%.0f)", median(sub$age, na.rm = TRUE),
                quantile(sub$age, .25, na.rm = TRUE), quantile(sub$age, .75, na.rm = TRUE)),
  male = sprintf("%.1f", 100 * mean(sub$gender.demographic == "male", na.rm = TRUE)),
  stIV = sprintf("%.1f", 100 * mean(sub$stage_iv, na.rm = TRUE)),
  tmb = sprintf("%.2f (%.2f-%.2f)", median(sub$TMB, na.rm = TRUE),
                quantile(sub$TMB, .25, na.rm = TRUE), quantile(sub$TMB, .75, na.rm = TRUE)),
  pdl1 = sprintf("%.2f", median(as.numeric(tc$X["CD274", sub$sample]), na.rm = TRUE)),
  events = sprintf("%d (%.1f)", sum(sub$OS, na.rm = TRUE), 100 * mean(sub$OS, na.rm = TRUE)),
  fu = sprintf("%.1f", median(sub$OS.time, na.rm = TRUE) / 30.44))
tab1 <- data.table(variable = c("Patients, n", "Age, years", "Male, %", "Stage IV, %",
                                "TMB, mut/Mb", "CD274 (PD-L1), log2(FPKM+1)", "OS events, n (%)",
                                "Median follow-up, months"),
  immune_activated = mk(cl[pheno == 1]), immune_cold = mk(cl[pheno == 0]))
fwrite(tab1, file.path(TDIR, "table1_tcga_baseline.csv"))
tab2 <- perf[, .(Signature = signature, Cohort = cohort, Cancer = cancer, Platform = platform,
  N = n, Responders = n_resp, `ORR %` = round(100 * orr, 1), Genes = n_gene_used,
  AUC = round(auc, 3), `95% CI` = sprintf("%.3f-%.3f", auc_lo, auc_hi),
  `OR/SD` = sprintf("%.2f (%.2f-%.2f)", or, or_lo, or_hi),
  `OS HR/SD` = ifelse(is.na(os_hr), "\u2014", sprintf("%.2f (P=%.3g)", os_hr, os_p)))]
fwrite(tab2, file.path(TDIR, "table2_cohort_performance.csv"))
tab3 <- pool_all[, .(Signature = signature, k, `Total N` = total_n,
  `Pooled AUC` = round(pooled_auc, 3), `95% CI` = sprintf("%.3f-%.3f", lo, hi),
  I2 = sprintf("%.1f%%", I2), `Q P` = signif(Qp, 3),
  `Pooled OR/SD` = sprintf("%.2f (%.2f-%.2f), P=%.3g", pooled_or, or_lo, or_hi, or_p))]
fwrite(tab3, file.path(TDIR, "table3_signature_comparison.csv"))
cat("\n===== Table 1 =====\n"); print(tab1)
cat("\n===== Table 3 =====\n"); print(tab3)

## ---------- 8. 汇总 RDS ----------
saveRDS(list(perf = perf, pool = pool_all, loco = loco, beta = beta,
             perm = list(auc_obs = auc_obs, p = p_perm, median = median(pa)),
             disc = disc, cluster_sizes = disc$cluster_sizes),
        file.path(RDIR, "final_results.rds"))
msg("FINAL DONE")
