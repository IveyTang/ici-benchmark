## ============================================================
## 08_jtm_tables.R
## 为 JTM 投稿版生成两张表的权威内容：
##   (1) Table 2  —— IAS 26 基因及其 LASSO 系数（从 final_results.rds 的 beta 提取）
##   (2) Table 1  —— TCGA-HNSC 基线特征 + 完整 P 值（不编造，全部实算）
## 所有数值均来自已冻结的分析对象，无手工输入。
## ============================================================
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })

fin <- readRDS(file.path(RDIR, "final_results.rds"))
tc  <- readRDS(file.path(RDIR, "tcga_hnsc.rds"))
msg("final_results.rds 元素：", paste(names(fin), collapse = ", "))

## ---------- (1) Table 2：26 基因 IAS 系数 ----------
bv <- fin$beta
msg("beta 长度 = ", length(bv), "  非零 = ", sum(bv != 0))
tb2 <- data.table(gene = names(bv), coefficient = as.numeric(bv))
tb2 <- tb2[coefficient != 0]
tb2[, abscoef := abs(coefficient)]
setorder(tb2, -abscoef)
tb2[, rank_abs := .I]
tb2[, coefficient := round(coefficient, 4)]
fwrite(tb2, file.path(TDIR, "jtm_table2_ias_genes.csv"))
msg("Table 2 已写出，基因数 = ", nrow(tb2))
print(tb2)

## ---------- (2) Table 1：基线特征 + P 值 ----------
cl <- as.data.table(tc$clin)
cl[, age   := as.numeric(age_at_index.demographic)]
cl[, stage_iv := grepl("IV", ajcc_pathologic_stage.diagnoses)]
## 表型标签来自 06_final.R 落盘的 tcga_ias_scores.csv（sample, IAS, phenotype, split）
ph <- fread(file.path(TDIR, "tcga_ias_scores.csv"))[, .(sample, pheno = phenotype)]
cl <- merge(cl, ph, by = "sample")
cl <- merge(cl, as.data.table(tc$tmb)[, .(sample, TMB)], by = "sample", all.x = TRUE)
cl <- merge(cl, as.data.table(tc$surv), by = "sample", all.x = TRUE)
cl[, cd274 := as.numeric(tc$X["CD274", sample])]
cl[, male := gender.demographic == "male"]

a <- cl[pheno == 1]; c0 <- cl[pheno == 0]
fm <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x, na.rm = TRUE),
                          quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
pk <- function(p) ifelse(p < 2.2e-16, "< 2.2e-16", formatC(p, format = "e", digits = 1))
fu <- function(x) sprintf("%.1f", median(x, na.rm = TRUE) / 30.44)

p_age  <- wilcox.test(a$age, c0$age)$p.value
p_sex  <- chisq.test(table(cl$pheno, cl$male))$p.value
p_stiv <- chisq.test(table(cl$pheno, cl$stage_iv))$p.value
p_tmb  <- wilcox.test(a$TMB, c0$TMB)$p.value
p_pdl1 <- wilcox.test(a$cd274, c0$cd274)$p.value
lung   <- survdiff(Surv(OS.time, OS) ~ pheno, data = cl[!is.na(OS.time) & !is.na(OS)])
p_os   <- 1 - pchisq(lung$chisq, length(lung$n) - 1)

tb1 <- data.table(
  variable = c("Patients, n", "Age, years, median (IQR)", "Male, %",
               "Stage IV, %", "TMB, mut/Mb, median (IQR)",
               "CD274 (PD-L1), log2(FPKM+1), median",
               "OS events, n (%)", "Median follow-up, months"),
  immune_activated = c(nrow(a), sprintf("%.0f (%.0f\u2013%.0f)", median(a$age, na.rm=TRUE),
                        quantile(a$age,.25,na.rm=TRUE), quantile(a$age,.75,na.rm=TRUE)),
                        sprintf("%.1f", 100*mean(a$male, na.rm=TRUE)),
                        sprintf("%.1f", 100*mean(a$stage_iv, na.rm=TRUE)),
                        fm(a$TMB), sprintf("%.2f", median(a$cd274, na.rm=TRUE)),
                        sprintf("%d (%.1f)", sum(a$OS, na.rm=TRUE), 100*mean(a$OS, na.rm=TRUE)),
                        fu(a$OS.time)),
  immune_cold = c(nrow(c0), sprintf("%.0f (%.0f\u2013%.0f)", median(c0$age, na.rm=TRUE),
                        quantile(c0$age,.25,na.rm=TRUE), quantile(c0$age,.75,na.rm=TRUE)),
                        sprintf("%.1f", 100*mean(c0$male, na.rm=TRUE)),
                        sprintf("%.1f", 100*mean(c0$stage_iv, na.rm=TRUE)),
                        fm(c0$TMB), sprintf("%.2f", median(c0$cd274, na.rm=TRUE)),
                        sprintf("%d (%.1f)", sum(c0$OS, na.rm=TRUE), 100*mean(c0$OS, na.rm=TRUE)),
                        fu(c0$OS.time)),
  P = c("\u2014", pk(p_age), pk(p_sex), pk(p_stiv), pk(p_tmb), pk(p_pdl1), pk(p_os), "\u2014"))
fwrite(tb1, file.path(TDIR, "jtm_table1_baseline_pvalues.csv"))
msg("Table 1 已写出（含实算 P 值）")
print(tb1)
msg("JTM TABLES DONE")
