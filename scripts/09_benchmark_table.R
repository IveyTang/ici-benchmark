## ============================================================
## 09_benchmark_table.R — 27 个打分构建策略的基准总表
##
## 输入：exp_ias_results.rds、exp_ias2_results.rds、
##       exp_ias3_results.rds、exp_ias4_results.rds
##       （由 exp_ias.R / exp_ias2.R / exp_ias3.R / exp_ias4.R 生成）
## 输出：results/tables/table4_benchmark_all_strategies.csv
##
## 每个策略的合并 AUC 均以同一随机效应模型（DerSimonian-Laird，
## logit 尺度 + Hanley-McNeil 方差）在相同的 6 个队列上计算，故各行可比。
## ============================================================
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))

e1 <- readRDS(file.path(RDIR, "exp_ias_results.rds"))$res      # 集合级方案
e2 <- readRDS(file.path(RDIR, "exp_ias2_results.rds"))$R       # 基因级方案
e3 <- readRDS(file.path(RDIR, "exp_ias3_results.rds"))$R       # IFN-gamma 锚点派生方案
e4 <- readRDS(file.path(RDIR, "exp_ias4_results.rds"))$S       # 跨队列留一训练方案

A <- rbindlist(list(
  e1[, .(strategy = variant, n_feature = n_sig,  pooled_auc, lo, hi, I2)],
  e2[, .(strategy = variant, n_feature = n_gene, pooled_auc, lo, hi, I2)],
  e3[, .(strategy = variant, n_feature = n_gene, pooled_auc, lo, hi, I2)],
  ## e4 只贡献参照签名与 IMvigor210 训练方案；LOCO_* 留一训练方案另表报告
  ## 注意 e4 中训练方案的原始名称是 "TrainIMvigor210"（07 脚本中另起别名 TRAIN_IMV）
  e4[grepl("^REF_|^Train", scheme),
     .(strategy = scheme, n_feature = n_gene, pooled_auc, lo, hi, I2)]))
A <- unique(A, by = "strategy")

## 与 07_figures_extra.R 的森林图保持一致：剔除被参照签名覆盖的冗余方案
##   F_IFNG_6  —— 与 REF_IFNG6 属同一签名，仅因数据源不同而重复出现
##   REF_CYT2  —— 2 基因细胞溶解活性签名，效应量估计不稳定
A <- A[!strategy %in% c("F_IFNG_6", "REF_CYT2", "REF_IFNG_6")]

## 层次归属采用纯规则，无人工赋值：
##   V*               -> 基因集富集分数级（GSVA / ssGSEA 派生）
##   TrainIMvigor210  -> 在外部队列上训练的方案
##   其余             -> 基因表达级（含已发表参照签名）
A[, level := fifelse(grepl("^V[0-9]", strategy), "Gene-set enrichment score",
              fifelse(strategy == "TrainIMvigor210", "Published reference / trained",
                      "Gene expression"))]
setorder(A, -pooled_auc)

## 统一输出精度与列序，便于逐字比对
A[, `:=`(pooled_auc = round(pooled_auc, 4), lo = round(lo, 4),
         hi = round(hi, 4), I2 = round(I2, 1))]
setcolorder(A, c("strategy", "level", "n_feature", "pooled_auc", "lo", "hi", "I2"))

fwrite(A, file.path(TDIR, "table4_benchmark_all_strategies.csv"))
msg("table4_benchmark_all_strategies.csv 完成 — ", nrow(A), " 个策略")
print(A[, .(strategy, level, n_feature, pooled_auc, lo, hi)])

## 两个层次的合并 AUC 区间（正文核心结论的直接依据）
msg("按层次汇总:")
print(A[, .(min = min(pooled_auc), max = max(pooled_auc), n = .N, max_I2 = max(I2)),
        by = level])
