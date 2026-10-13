#!/usr/bin/env bash
# ncbi_geo.sh — 绕过域名级封锁访问 NCBI/GEO
# 两个可用通道：
#   1) FTP  : ftp://ftp.ncbi.nlm.nih.gov/...            （最简，优先）
#   2) HTTPS: --resolve <host>:443:<NCBI_IP> --noproxy '*' -k  （回退）
# 用法:
#   bash ncbi_geo.sh ls   gse159067            # 列出 suppl 目录
#   bash ncbi_geo.sh get  gse159067/suppl/xxx.gz  out.gz
set -u
NCBI_IP="${NCBI_IP:-130.14.250.31}"
HOST="ftp.ncbi.nlm.nih.gov"
BASE_PATH="geo/series"

acc_to_dir() {              # GSE159067 -> GSE159nnn
  local a; a=$(echo "$1" | tr '[:lower:]' '[:upper:]')
  local digits=${a#GSE}
  local n=${#digits}
  local prefix=${digits:0:$((n-3))}
  echo "GSE${prefix}nnn/$a"
}

cmd="${1:-}"
case "$cmd" in
  ls)
    d=$(acc_to_dir "$2")
    url="ftp://$HOST/$BASE_PATH/$d/suppl/"
    curl -s --connect-timeout 20 -m 60 "$url" 2>/dev/null | sed 's/.*href="//;s/".*//' | grep -v '^/' || {
      # FTP 列目录失败则用 HTTPS+resolve
      curl -s -k --noproxy '*' --resolve "$HOST:443:$NCBI_IP" -m 60 "https://$HOST/$BASE_PATH/$d/suppl/" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//' | grep -v '^/'
    }
    ;;
  lsmat)
    d=$(acc_to_dir "$2")
    curl -s -k --noproxy '*' --resolve "$HOST:443:$NCBI_IP" -m 60 "https://$HOST/$BASE_PATH/$d/matrix/" 2>/dev/null | grep -oE 'href="GSE[^"]+"' | sed 's/href="//;s/"//'
    ;;
  get)
    rel="$2"; out="$3"
    mkdir -p "$(dirname "$out")"
    # 先试 FTP
    if curl -s --connect-timeout 20 -m 1800 -o "$out" "ftp://$HOST/$BASE_PATH/$rel" 2>/dev/null && [ -s "$out" ]; then
      echo "OK(ftp) $out $(stat -c%s "$out" 2>/dev/null || wc -c <"$out") bytes"; exit 0
    fi
    # 回退 HTTPS+resolve
    if curl -s -k --noproxy '*' --resolve "$HOST:443:$NCBI_IP" -m 1800 -o "$out" "https://$HOST/$BASE_PATH/$rel" 2>/dev/null && [ -s "$out" ]; then
      echo "OK(https) $out $(stat -c%s "$out" 2>/dev/null || wc -c <"$out") bytes"; exit 0
    fi
    echo "FAIL $rel"; exit 1
    ;;
  *)
    echo "用法: bash ncbi_geo.sh {ls|lsmat|get} ..."; exit 1;;
esac
