## ============================================================
## 00_config.R — shared configuration
##
## 仓库根目录 = 含 data/ 与 results/ 的目录。
## 解析顺序：
##   1) 环境变量 ICI_WORK（推荐显式指定）
##   2) 当前工作目录
## 因此从仓库根目录运行 scripts/*.R 即可，无需改任何路径。
##
## 目录约定
##   data/            原始数据（需先运行 scripts/download_data.sh）
##   results/tables/  派生表（CSV）
##   results/figures/ 派生图（PNG）
##   results/rds/     中间对象（RDS，体积大者不纳入版本控制）
## ============================================================

WORK <- Sys.getenv("ICI_WORK", unset = getwd())
if (!dir.exists(WORK)) stop("ICI_WORK 指向的目录不存在: ", WORK)

DL   <- file.path(WORK, "data")
OUT  <- file.path(WORK, "results")
RDIR <- file.path(OUT, "rds")
TDIR <- file.path(OUT, "tables")
FDIR <- file.path(OUT, "figures")
for (d in c(OUT, RDIR, TDIR, FDIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## 可选：本地包库。优先环境变量 ICI_LIB，其次 <root>/Rlib。
## 两者都不存在时使用标准 .libPaths()（即正常安装在系统库/用户库中的包）。
.lib_cand <- c(Sys.getenv("ICI_LIB", unset = ""), file.path(WORK, "Rlib"))
.lib_cand <- .lib_cand[nzchar(.lib_cand) & dir.exists(.lib_cand)]
if (length(.lib_cand)) .libPaths(c(.lib_cand, .libPaths()))

options(timeout = 3600, stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

## 统一的基因集富集打分。
## 注意：发现阶段（02_discovery.R）显式指定 method = "gsva"。
## ssGSEA 仅用于敏感性分析；其样本内秩基性质会抹掉免疫转录本总丰度信号，
## 使样本间 Spearman 相关中位数高达约 0.93，并令一致性聚类退化为 517/3。
run_enrich <- function(X, gsc, method = "ssgsea") {
  if (method == "ssgsea") {
    r <- tryCatch(suppressWarnings(GSVA::gsva(GSVA::ssgseaParam(X, gsc, normalize = TRUE))),
                  error = function(e) NULL)
    if (!is.null(r)) return(r)
    r <- tryCatch(suppressWarnings(GSVA::gsva(X, gsc, method = "ssgsea",
                                              kcdf = "Gaussian", verbose = FALSE)),
                  error = function(e) NULL)
    if (!is.null(r)) return(r)
  }
  r <- tryCatch(suppressWarnings(GSVA::gsva(GSVA::gsvaParam(X, gsc, kcdf = "Gaussian",
                                                            maxDiff = TRUE, verbose = FALSE))),
                error = function(e) NULL)
  if (!is.null(r)) return(r)
  suppressWarnings(GSVA::gsva(X, gsc, method = "gsva", kcdf = "Gaussian",
                              mx.diff = TRUE, verbose = FALSE))
}
