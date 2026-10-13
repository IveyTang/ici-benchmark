#!/usr/bin/env bash
# ============================================================
# run_all.sh — 按顺序执行完整分析管线
#
# 用法：
#   bash scripts/run_all.sh            # 从 01_prep 跑到 08
#   bash scripts/run_all.sh 06         # 从 06_final 开始（复用已有缓存）
#
# 前置：数据已就位（bash scripts/download_data.sh）
# 日志：写入 logs/，同时输出到终端
# ============================================================
set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export ICI_WORK="${ICI_WORK:-$ROOT}"
mkdir -p logs

RSCRIPT="${RSCRIPT:-Rscript}"
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "找不到 Rscript，请设置 RSCRIPT=/path/to/Rscript"; exit 1; }

# 执行顺序：准备 -> 发现 -> 四组探索 -> 标记面板 -> 最终验证 -> 图 -> 基准表 -> 稿件用表
CHAIN=(
  "01_prep.R"           # 队列与表达矩阵
  "02_discovery.R"      # TCGA-HNSC 免疫表型 + 派生分数
  "exp_ias.R"           # 探索 1：集合级方案
  "exp_ias2.R"          # 探索 2：基因级方案
  "exp_ias3.R"          # 探索 3：IFN-gamma 锚点派生方案
  "exp_ias4.R"          # 探索 4：跨队列留一训练
  "04_markers.R"        # 免疫细胞标记面板
  "06_final.R"          # 跨队列验证（主结果）
  "07_figures_extra.R"  # fig6 / fig7
  "09_benchmark_table.R" # 27 方案基准表
  "08_jtm_tables.R"     # 稿件用表
)

START="${1:-}"
RUN=1
[ -n "$START" ] && RUN=0

for s in "${CHAIN[@]}"; do
  [ -z "$START" ] || case "$s" in "$START"*) RUN=1;; esac
  [ "$RUN" -eq 1 ] || continue
  echo "=============================================================="
  echo ">> $s    $(date '+%Y-%m-%d %H:%M:%S')"
  echo "=============================================================="
  if ! "$RSCRIPT" "scripts/$s" 2>&1 | tee "logs/${s%.R}.log"; then
    echo "!! $s 失败，管线中止"
    exit 1
  fi
done

echo
echo "全部完成。结果位于 results/（tables/ 表，figures/ 图，rds/ 中间对象）"
