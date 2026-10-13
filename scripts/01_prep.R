## ============================================================
## 01_prep.R — 数据准备：TCGA-HNSC 发现队列 + 6 个 ICI 验证队列
## 含缓存：重复运行会跳过耗时的注释解析与 TCGA 读取
## ============================================================
## ---- 可移植前导：定位本脚本所在目录并加载共享配置 ----
.local_this <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.local_this)) dirname(normalizePath(.local_this[1])) else getwd()
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(openxlsx) })

## ---------- A. 基因注释（带缓存）----------
f_ens <- file.path(RDIR, "map_ensembl_symbol.rds")
f_ent <- file.path(RDIR, "map_entrez_symbol.rds")
if (file.exists(f_ens) && file.exists(f_ent)) {
  msg("加载缓存的基因注释")
  ens_map <- readRDS(f_ens); ent_map <- readRDS(f_ent)
} else {
  msg("解析 NCBI gene_info 注释")
  gi <- fread(file.path(DL, "anno/Homo_sapiens.gene_info.gz"), sep = "\t",
              header = TRUE, quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gi, 1, "tax_id")
  gi <- gi[tax_id == 9606]
  ens_map <- unique(rbindlist(lapply(seq_len(nrow(gi)), function(i) {
    m <- regmatches(gi$dbXrefs[i], gregexpr("Ensembl:ENSG[0-9]+", gi$dbXrefs[i], perl = TRUE))[[1]]
    if (!length(m)) return(NULL)
    data.table(ens = sub("^Ensembl:", "", m), symbol = gi$Symbol[i])
  })), by = "ens")
  ent_map <- unique(data.table(entrez = as.character(gi$GeneID), symbol = gi$Symbol), by = "entrez")
  saveRDS(ens_map, f_ens); saveRDS(ent_map, f_ent)
}
msg("  Ensembl->Symbol: ", nrow(ens_map), "  Entrez->Symbol: ", nrow(ent_map))

## ---------- B. TCGA-HNSC（带缓存）----------
f_tcga <- file.path(RDIR, "tcga_hnsc.rds")
if (file.exists(f_tcga)) {
  msg("加载缓存的 TCGA-HNSC")
  tc <- readRDS(f_tcga); Xtcga <- tc$X; surv <- as.data.table(tc$surv)
  clin <- as.data.table(tc$clin); tmb <- as.data.table(tc$tmb)
} else {
  msg("读取 TCGA-HNSC 表达矩阵 (Xena GDC star_fpkm)")
  tcga <- fread(file.path(DL, "xena/TCGA-HNSC.star_fpkm.tsv.gz"), sep = "\t",
                header = TRUE, showProgress = FALSE)
  setnames(tcga, 1, "ens")
  tcga[, ens := sub("\\..*$", "", ens)]
  tcga <- merge(ens_map[, .(ens, symbol)], tcga, by = "ens")
  tcga[, ens := NULL]
  ## 快得多的聚合：先转矩阵再按 symbol 分组求均值
  smps <- setdiff(names(tcga), "symbol")
  M <- as.matrix(tcga[, ..smps]); storage.mode(M) <- "double"
  grp <- tcga$symbol
  Mc <- rowsum(M, group = grp, reorder = TRUE)
  n <- as.vector(rowsum(matrix(1, nrow(M), 1), group = grp, reorder = TRUE))
  Ma <- Mc / n
  rm(M, Mc, tcga); invisible(gc())
  keep <- substr(colnames(Ma), 14, 15) == "01"      # 仅原发肿瘤
  Xtcga <- Ma[, keep, drop = FALSE]
  Xtcga <- Xtcga[rowSums(!is.finite(Xtcga)) == 0, , drop = FALSE]
  vv <- rowVars(Xtcga); Xtcga <- Xtcga[vv > 0, , drop = FALSE]
  msg("  TCGA expr: ", nrow(Xtcga), " genes x ", ncol(Xtcga), " primary tumours")

  surv <- fread(file.path(DL, "xena/TCGA-HNSC.survival.tsv.gz"), sep = "\t", header = TRUE)
  setnames(surv, c("sample", "OS.time", "OS", "patient"))
  surv <- surv[sample %in% colnames(Xtcga)]
  surv[, OS.time := as.numeric(OS.time)]; surv[, OS := as.integer(OS)]

  clin <- fread(file.path(DL, "xena/TCGA-HNSC.clinical.tsv.gz"), sep = "\t",
                header = TRUE, fill = TRUE, showProgress = FALSE)
  clin <- clin[sample %in% colnames(Xtcga)]
  keepc <- intersect(c("sample", "gender.demographic", "age_at_index.demographic",
                       "ajcc_pathologic_stage.diagnoses", "primary_site",
                       "vital_status.demographic"), names(clin))
  clin <- clin[, keepc, with = FALSE]

  mut <- fread(file.path(DL, "xena/TCGA-HNSC.somaticmutation_wxs.tsv.gz"), sep = "\t",
               header = TRUE, showProgress = FALSE)
  tmb <- mut[, .(n_mut = .N), by = sample][sample %in% colnames(Xtcga)]
  tmb[, TMB := n_mut / 38]
  saveRDS(list(X = Xtcga, surv = surv, clin = clin, tmb = tmb), f_tcga)
}
msg("  生存: ", nrow(surv), " 例, 事件 ", sum(surv$OS), ", 中位随访 ",
    round(median(surv$OS.time) / 30.44, 1), " 月; TMB 中位 ",
    round(median(tmb$TMB), 2), " mut/Mb")

## ---------- C. 工具函数 ----------
## 解析 series matrix 样本元数据。
## 关键：GEO 的 characteristics 行会因样本特征数不一致而错位（同一行内混入别的标签），
## 因此必须逐单元格解析，让每个值自带标签，再按 (id, label) 透视。
series_matrix_meta <- function(gzfile, id_field = "Sample_description") {
  ln <- readLines(gzfile(gzfile), warn = FALSE)
  ln <- ln[grepl("^!(Sample_title|Sample_geo_accession|Sample_description|Sample_characteristics)", ln)]
  vals_of <- function(prefix) {
    l <- ln[startsWith(ln, prefix)][1]
    if (is.na(l)) return(NULL)
    gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  }
  ids <- vals_of(paste0("!", id_field, "\t"))
  if (is.null(ids)) ids <- vals_of("!Sample_title\t")
  ids <- trimws(ids)
  n <- length(ids)
  sid <- character(0); lab <- character(0); vl <- character(0)
  for (l in ln[startsWith(ln, "!Sample_characteristics")]) {
    v <- gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
    for (i in seq_along(v)) {
      cell <- v[i]
      if (is.na(cell) || !nzchar(cell) || i > n || !grepl(":", cell, fixed = TRUE)) next
      sid <- c(sid, ids[i])
      lab <- c(lab, gsub("[^A-Za-z0-9]+", ".", trimws(sub(":.*$", "", cell))))
      vl  <- c(vl, sub("^[^:]*:\\s*", "", cell))
    }
  }
  d <- unique(data.table(id = sid, label = lab, value = vl), by = c("id", "label"))
  wide <- dcast(d, id ~ label, value.var = "value")
  ## 保证样本顺序与 series matrix 一致
  wide <- wide[match(ids, id)]
  wide
}
## 用正则容错取字段（返回与 ids 对齐的向量；不存在则 NULL）
mf <- function(meta, pattern, ids = NULL) {
  hit <- grep(pattern, names(meta), ignore.case = TRUE, value = TRUE)
  if (!length(hit)) { warning("字段未找到: ", pattern); return(NULL) }
  v <- meta[[hit[1]]]
  if (is.null(ids)) v else v[match(ids, meta$id)]
}

## ---------- D. 各 ICI 队列 ----------
cohorts <- list()

## D1. GSE159067 — HNSCC, HTG EdgeSeq (n=102)【HNSCC 专属验证队列】
msg("读取 GSE159067 (HNSCC, HTG EdgeSeq)")
m <- fread(file.path(DL, "geo/GSE159067_IHN_log2cpm_data.txt.gz"), sep = "\t",
           header = TRUE, showProgress = FALSE)
setnames(m, 1, "gene")
m <- m[!grepl("^(NEG|POS)_CTRL|^(NEG|POS)_", gene)]
X <- as.matrix(m[, -1]); rownames(X) <- m$gene; storage.mode(X) <- "double"
md <- series_matrix_meta(file.path(DL, "geo/GSE159067_series_matrix.txt.gz"), "Sample_description")
ids <- colnames(X); sel <- match(ids, md$id)
meta <- data.table(
  id = ids, site = mf(md, "^site")[sel], age = as.numeric(mf(md, "^age$")[sel]),
  sex = mf(md, "^Sex$")[sel], hpv = mf(md, "^hpv")[sel],
  line = mf(md, "^immunotherapy\\.line")[sel],
  best = mf(md, "^best\\.response")[sel],
  os_time = as.numeric(mf(md, "^os\\.month")[sel]),
  os_event = as.integer(mf(md, "^os\\.event")[sel]),
  pfs_time = as.numeric(mf(md, "^pfs\\.month")[sel]),
  pfs_event = as.integer(mf(md, "^pfs\\.event")[sel]))
meta[, response := fifelse(grepl("^(CR|PR)$", best), 1L,
                    fifelse(grepl("^(SD|PD)$", best), 0L, NA_integer_))]
cohorts[["GSE159067"]] <- list(name = "GSE159067", cancer = "HNSCC",
  platform = "HTG EdgeSeq (targeted)", X = X, meta = meta, response = meta$response,
  os_time = meta$os_time, os_event = meta$os_event, pfs_time = meta$pfs_time, pfs_event = meta$pfs_event)
msg("  样本 ", ncol(X), "，应答 ", sum(meta$response == 1, na.rm = TRUE),
    " / 非应答 ", sum(meta$response == 0, na.rm = TRUE))
print(table(meta$best, useNA = "ifany"))

## D2. GSE78220 — melanoma, anti-PD-1 (Hugo 2016, n=28)
msg("读取 GSE78220 (melanoma)")
fp <- as.data.table(read.xlsx(file.path(DL, "geo/GSE78220_PatientFPKM.xlsx"), sheet = "FPKM"))
setnames(fp, 1, "gene")
X <- as.matrix(fp[, -1]); rownames(X) <- fp$gene; storage.mode(X) <- "double"
colnames(X) <- sub("\\.baseline$", "", colnames(X))
X <- log2(X + 1)
md <- series_matrix_meta(file.path(DL, "geo/GSE78220_series_matrix.txt.gz"), "Sample_title")
sel <- match(colnames(X), md$id)
meta <- data.table(id = colnames(X), resp_raw = mf(md, "^anti\\.pd\\.1\\.response")[sel])
meta[, response := fifelse(grepl("Complete Response|Partial Response", resp_raw), 1L,
                    fifelse(grepl("Progressive Disease|Stable Disease", resp_raw), 0L, NA_integer_))]
cohorts[["GSE78220"]] <- list(name = "GSE78220", cancer = "Melanoma",
  platform = "RNA-seq (FPKM)", X = X, meta = meta, response = meta$response)
msg("  样本 ", ncol(X), "，应答 ", sum(meta$response == 1, na.rm = TRUE),
    " / 非应答 ", sum(meta$response == 0, na.rm = TRUE))

## D3. GSE91061 — melanoma, nivolumab (Riaz 2017, 仅基线样本)
msg("读取 GSE91061 (melanoma)")
m <- fread(file.path(DL, "geo/GSE91061_fpkm.csv.gz"), sep = ",", header = TRUE, showProgress = FALSE)
setnames(m, 1, "entrez"); m <- m[!is.na(entrez) & as.character(entrez) != ""]
m[, entrez := as.character(entrez)]
m <- merge(ent_map, m, by = "entrez"); m[, entrez := NULL]
sm <- setdiff(names(m), "symbol")
M <- as.matrix(m[, ..sm]); storage.mode(M) <- "double"
Mc <- rowsum(M, group = m$symbol, reorder = TRUE)
nn <- as.vector(rowsum(matrix(1, nrow(M), 1), group = m$symbol, reorder = TRUE))
X <- log2(Mc / nn + 1); rm(m, M, Mc); invisible(gc())
md <- series_matrix_meta(file.path(DL, "geo/GSE91061_series_matrix.txt.gz"), "Sample_title")
sel <- match(colnames(X), md$id)
meta <- data.table(id = colnames(X), visit = mf(md, "^visit")[sel],
                   resp_raw = mf(md, "^response$")[sel])
meta[, response := fifelse(grepl("^PRCR$", resp_raw), 1L,
                    fifelse(grepl("^(PD|SD)$", resp_raw), 0L, NA_integer_))]
pre <- !is.na(meta$visit) & meta$visit == "Pre"
X <- X[, pre, drop = FALSE]; meta <- meta[pre]
cohorts[["GSE91061"]] <- list(name = "GSE91061", cancer = "Melanoma",
  platform = "RNA-seq (FPKM)", X = X, meta = meta, response = meta$response)
msg("  样本 ", ncol(X), "（仅基线），应答 ", sum(meta$response == 1, na.rm = TRUE),
    " / 非应答 ", sum(meta$response == 0, na.rm = TRUE))

## D4. GSE135222 — NSCLC, anti-PD-1/PD-L1 (Jung 2019, n=27)
msg("读取 GSE135222 (NSCLC)")
m <- fread(file.path(DL, "geo/GSE135222_exp.tsv.gz"), sep = "\t", header = TRUE, showProgress = FALSE)
setnames(m, 1, "ens"); m[, ens := sub("\\..*$", "", ens)]
m <- merge(ens_map[, .(ens, symbol)], m, by = "ens"); m[, ens := NULL]
sm <- setdiff(names(m), "symbol")
M <- as.matrix(m[, ..sm]); storage.mode(M) <- "double"
Mc <- rowsum(M, group = m$symbol, reorder = TRUE)
nn <- as.vector(rowsum(matrix(1, nrow(M), 1), group = m$symbol, reorder = TRUE))
X <- log2(Mc / nn + 1); rm(m, M, Mc); invisible(gc())
md <- series_matrix_meta(file.path(DL, "geo/GSE135222_series_matrix.txt.gz"), "Sample_title")
md <- copy(md); md[, idn := gsub(" ", "", id)]
sel <- match(gsub(" ", "", colnames(X)), md$idn)
meta <- data.table(id = colnames(X),
  pfs_time = as.numeric(mf(md, "^pfs\\.time")[sel]),
  pfs_event = as.integer(mf(md, "^progression\\.free")[sel]))
## pfs.time 单位为「天」（最大 618 天 ≈ 20 个月）；持久获益定义为 PFS >= 6 个月
meta[, response := fifelse(!is.na(pfs_time) & pfs_time >= 182, 1L,
                    fifelse(!is.na(pfs_time), 0L, NA_integer_))]
cohorts[["GSE135222"]] <- list(name = "GSE135222", cancer = "NSCLC",
  platform = "RNA-seq (counts)", X = X, meta = meta, response = meta$response,
  pfs_time = meta$pfs_time, pfs_event = meta$pfs_event)
msg("  样本 ", ncol(X), "，获益(PFS>=6m) ", sum(meta$response == 1), " / 非获益 ",
    sum(meta$response == 0))

## D5. GSE176307 — 尿路上皮癌 (BACI), anti-PD-1
msg("读取 GSE176307 (bladder/BACI)")
m <- fread(file.path(DL, "geo/GSE176307_norm.csv.gz"), sep = ",", header = TRUE, showProgress = FALSE)
setnames(m, 1, "sample"); gn <- setdiff(names(m), "sample")
X <- t(as.matrix(m[, ..gn])); storage.mode(X) <- "double"; colnames(X) <- m$sample
md <- series_matrix_meta(file.path(DL, "geo/GSE176307_series_matrix.txt.gz"), "Sample_title")
md <- copy(md); md[, baci := sub(".*(BACI[0-9]+).*", "\\1", id)]
sel <- match(colnames(X), md$baci)
meta <- data.table(id = colnames(X), resp_raw = mf(md, "^io\\.response")[sel],
  pfs_time = as.numeric(mf(md, "^pfs$")[sel]),
  progressed = as.integer(mf(md, "^progressed")[sel]),
  tmb = as.numeric(mf(md, "^tmb$")[sel]))
meta[, response := fifelse(grepl("CR|PR", resp_raw, ignore.case = TRUE), 1L,
                    fifelse(grepl("PD|SD", resp_raw, ignore.case = TRUE), 0L, NA_integer_))]
cohorts[["GSE176307"]] <- list(name = "GSE176307", cancer = "Urothelial (BACI)",
  platform = "RNA-seq (log-normalised)", X = X, meta = meta, response = meta$response,
  tmb = meta$tmb, pfs_time = meta$pfs_time, pfs_event = meta$progressed)
msg("  样本 ", ncol(X), "，应答 ", sum(meta$response == 1, na.rm = TRUE),
    " / 非应答 ", sum(meta$response == 0, na.rm = TRUE), " / 缺失 ", sum(is.na(meta$response)))
print(table(meta$resp_raw, useNA = "ifany"))

## D6. IMvigor210 — 尿路上皮癌, atezolizumab (Mariathasan 2018, n=348 / 298 可评估)
## cds 为 DESeq v1 的 CountDataSet；DESeq 未安装，故直接读取 S4 属性槽位取原始 counts
f_imv <- file.path(RDIR, "imvigor210_raw.rds")
if (file.exists(f_imv)) {
  msg("加载缓存的 IMvigor210")
  im <- readRDS(f_imv); Xi <- im$X; pdi <- im$pd
} else {
  msg("读取 IMvigor210 (cds.RData, DESeq v1 CountDataSet)")
  e <- new.env()
  load(file.path(DL, "gh/IMvigor210/data/cds.RData"), envir = e)
  cds <- get("cds", envir = e)
  M <- get("counts", envir = attr(cds, "assayData"))          # Entrez x 样本, 原始 counts
  pdi <- as.data.frame(attr(attr(cds, "phenoData"), "data"))
  msg("  原始矩阵: ", nrow(M), " 基因 x ", ncol(M), " 样本")
  idx <- match(ent_map$entrez, rownames(M)); ok <- !is.na(idx)
  Mi <- M[idx[ok], , drop = FALSE]; sy <- ent_map$symbol[ok]
  Mc <- rowsum(Mi, group = sy, reorder = TRUE)
  nn <- as.vector(rowsum(matrix(1, nrow(Mi), 1), group = sy, reorder = TRUE))
  lib <- colSums(Mc); lib[lib <= 0] <- 1
  Xi <- log2(sweep(Mc, 2, lib, "/") * 1e6 + 1)               # log2 CPM
  colnames(Xi) <- colnames(M); rm(M, Mi, Mc, e, cds); invisible(gc())
  saveRDS(list(X = Xi, pd = pdi), f_imv)
}
br <- as.character(pdi$binaryResponse)
resp_i <- fifelse(br == "CR/PR", 1L, fifelse(br == "SD/PD", 0L, NA_integer_))
meta_i <- data.table(id = colnames(Xi), response = resp_i,
  os_time = as.numeric(pdi$os) * 30.4375,                     # 月 -> 天
  os_event = as.integer(pdi$censOS),                          # 1 = 死亡事件
  tmb = as.numeric(pdi[["FMOne mutation burden per MB"]]),
  immune_phenotype = as.character(pdi[["Immune phenotype"]]),
  ic_level = as.character(pdi[["IC Level"]]),
  tc_level = as.character(pdi[["TC Level"]]),
  tcga_subtype = as.character(pdi[["TCGA Subtype"]]))
cohorts[["IMvigor210"]] <- list(name = "IMvigor210", cancer = "Urothelial",
  platform = "RNA-seq (targeted panel, log2 CPM)", X = Xi, meta = meta_i,
  response = resp_i, os_time = meta_i$os_time, os_event = meta_i$os_event,
  tmb = meta_i$tmb)
msg("  样本 ", ncol(Xi), "，可评估 ", sum(!is.na(resp_i)), "，应答 ", sum(resp_i == 1, na.rm = TRUE),
    " / 非应答 ", sum(resp_i == 0, na.rm = TRUE),
    "；事件 ", sum(meta_i$os_event == 1, na.rm = TRUE))

saveRDS(cohorts, file.path(RDIR, "ici_cohorts.rds"))
msg("完成: ", length(cohorts), " 个 ICI 队列")
for (nm in names(cohorts))
  msg("  - ", nm, ": ", ncol(cohorts[[nm]]$X), " samples x ", nrow(cohorts[[nm]]$X), " genes")
msg("PREP DONE")
