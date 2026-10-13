#!/usr/bin/env bash
# ============================================================
# download_data.sh — 抓取本研究全部公开输入数据到 data/
#
# 用法：
#   bash scripts/download_data.sh            # 全部
#   bash scripts/download_data.sh xena       # 仅某一步：xena|msigdb|anno|geo|imvigor210
#
# 说明：
#   * 数据总量约 420 MB，其中 TCGA-HNSC 表达矩阵约 140 MB。
#   * 已存在的文件默认跳过；加 --force 重新下载。
#   * 若所在网络封锁 NCBI，GEO 部分会自动走 ncbi_geo.sh 的 FTP / IP 直连通道。
# ============================================================
set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DL="$ROOT/data"
FORCE=0
for a in "$@"; do [ "$a" = "--force" ] && FORCE=1; done
STEP="${1:-all}"
[ "$STEP" = "--force" ] && STEP="all"

STEPS=(xena msigdb anno geo imvigor210)
[ "$STEP" != "all" ] && STEPS=("$STEP")

XENA="https://gdc-hub.s3.us-east-1.amazonaws.com/download"
MSIG="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs"
NCBI_GENE="https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
IMV="http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz"

FAIL=0
get() {                    # get <url> <outfile>
  local url="$1" out="$2"
  if [ -s "$out" ] && [ "$FORCE" -eq 0 ]; then
    echo "  [skip] $(basename "$out")"; return 0
  fi
  mkdir -p "$(dirname "$out")"
  echo "  [get ] $(basename "$out")"
  if curl -fsSL --connect-timeout 30 -m 3600 -o "$out.part" "$url"; then
    mv -f "$out.part" "$out"; return 0
  fi
  rm -f "$out.part"; echo "  [FAIL] $url"; FAIL=$((FAIL+1)); return 1
}

for s in "${STEPS[@]}"; do

# ---------- TCGA-HNSC（UCSC Xena / GDC hub）----------
if [ "$s" = "xena" ]; then
  echo "== TCGA-HNSC（Xena GDC hub）=="
  for f in TCGA-HNSC.star_fpkm.tsv.gz \
           TCGA-HNSC.survival.tsv.gz \
           TCGA-HNSC.clinical.tsv.gz \
           TCGA-HNSC.somaticmutation_wxs.tsv.gz; do
    get "$XENA/$f" "$DL/xena/$f"
  done
fi

# ---------- MSigDB（Hallmark + C7）----------
if [ "$s" = "msigdb" ]; then
  echo "== MSigDB 2024.1.Hs =="
  get "$MSIG/h.all.v2024.1.Hs.symbols.gmt" "$DL/msigdb/h.all.v2024.1.Hs.symbols.gmt"
  get "$MSIG/c7.all.v2024.1.Hs.symbols.gmt" "$DL/msigdb/c7.all.v2024.1.Hs.symbols.gmt"
fi

# ---------- NCBI 基因注释 ----------
if [ "$s" = "anno" ]; then
  echo "== NCBI Homo_sapiens.gene_info =="
  if ! get "$NCBI_GENE" "$DL/anno/Homo_sapiens.gene_info.gz"; then
    echo "  -> 回退：经代理中转"
    if command -v python >/dev/null 2>&1 && [ -f "$ROOT/scripts/proxy_fetch.py" ]; then
      python "$ROOT/scripts/proxy_fetch.py" "$NCBI_GENE" "$DL/anno/Homo_sapiens.gene_info.gz" || true
    fi
  fi
fi

# ---------- GEO 五个队列 ----------
if [ "$s" = "geo" ]; then
  echo "== GEO 验证队列 =="
  GEO="$ROOT/scripts/ncbi_geo.sh"
  # 参数：<目录前缀>/<accession>/<子目录>/<文件名>  <本地路径>
  fetch_geo() {   # fetch_geo <acc> <subdir> <filename> <localpath>
    local acc="$1" sub="$2" fn="$3" lp="$4"
    [ -s "$lp" ] && [ "$FORCE" -eq 0 ] && { echo "  [skip] $(basename "$lp")"; return 0; }
    mkdir -p "$(dirname "$lp")"
    echo "  [get ] $fn"
    bash "$GEO" get "${acc%???}nnn/$acc/$sub/$fn" "$lp" || FAIL=$((FAIL+1))
  }

  # GSE159067 — HNSCC, HTG EdgeSeq
  fetch_geo GSE159067 matrix GSE159067_series_matrix.txt.gz "$DL/geo/GSE159067_series_matrix.txt.gz"
  fetch_geo GSE159067 suppl  GSE159067_IHN_log2cpm_data.txt.gz "$DL/geo/GSE159067_IHN_log2cpm_data.txt.gz"

  # GSE176307 — urothelial (BACI)
  fetch_geo GSE176307 matrix GSE176307_series_matrix.txt.gz "$DL/geo/GSE176307_series_matrix.txt.gz"
  fetch_geo GSE176307 suppl  GSE176307_norm.csv.gz          "$DL/geo/GSE176307_norm.csv.gz"

  # GSE91061 — melanoma (Riaz)
  fetch_geo GSE91061 matrix GSE91061_series_matrix.txt.gz "$DL/geo/GSE91061_series_matrix.txt.gz"
  fetch_geo GSE91061 suppl  GSE91061_fpkm.csv.gz          "$DL/geo/GSE91061_fpkm.csv.gz"

  # GSE78220 — melanoma (Hugo)
  fetch_geo GSE78220 matrix GSE78220_series_matrix.txt.gz "$DL/geo/GSE78220_series_matrix.txt.gz"
  fetch_geo GSE78220 suppl  GSE78220_PatientFPKM.xlsx     "$DL/geo/GSE78220_PatientFPKM.xlsx"

  # GSE135222 — NSCLC (Jung)
  fetch_geo GSE135222 matrix GSE135222_series_matrix.txt.gz "$DL/geo/GSE135222_series_matrix.txt.gz"
  fetch_geo GSE135222 suppl  GSE135222_exp.tsv.gz          "$DL/geo/GSE135222_exp.tsv.gz"
fi

# ---------- IMvigor210 ----------
if [ "$s" = "imvigor210" ]; then
  echo "== IMvigor210（IMvigor210CoreBiologies 1.0.0, CC-BY 3.0）=="
  TAR="$DL/gh/IMvigor210CoreBiologies_1.0.0.tar.gz"
  get "$IMV" "$TAR"
  if [ -s "$TAR" ]; then
    echo "  [untar] 提取 data/*.RData 与 datalist"
    TMP="$DL/gh/_extract"
    mkdir -p "$TMP" "$DL/gh/IMvigor210/data"
    tar -xzf "$TAR" -C "$TMP"
    for f in cds.RData dat19.RData dat25.RData dat57.RData fmone.RData color_palettes.RData datalist; do
      find "$TMP" -name "$f" -exec cp -f {} "$DL/gh/IMvigor210/data/" \; 2>/dev/null || true
    done
    rm -rf "$TMP"
    ls -la "$DL/gh/IMvigor210/data/" | head
    echo "  注：cds.RData 为 DESeq v1 CountDataSet；01_prep.R 直读其 S4 槽位，无需安装 DESeq 包。"
  fi
fi

done

echo
if [ "$FAIL" -gt 0 ]; then
  echo "完成，但有 $FAIL 个文件下载失败。请检查网络后重跑对应步骤：bash scripts/download_data.sh <step>"
  exit 1
fi
echo "全部完成。数据位于 $DL"
echo "下一步：Rscript scripts/01_prep.R"
