#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
proxy_fetch.py — 通过公共 CORS 代理中转访问被网络封锁的资源（NCBI/GEO/zenodo 等）。

用法:
    python proxy_fetch.py <url> <outfile> [--text]
    python proxy_fetch.py --list <url>            # 只打印内容（文本）
"""
import sys, os, time, urllib.parse, urllib.request, gzip, io

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")

# 代理模板：{u} = URL 编码后的目标地址，{r} = 原始目标地址
PROXIES = [
    ("allorigins-raw",  "https://api.allorigins.win/raw?url={u}"),
    ("codetabs",        "https://api.codetabs.com/v1/proxy?quest={u}"),
    ("corsproxy.io",    "https://corsproxy.io/?url={u}"),
    ("thingproxy",      "https://thingproxy.freeboard.io/fetch/{r}"),
    ("allorigins-get",  "https://api.allorigins.win/get?url={u}"),
]


ONLY = os.environ.get("PROXY_ONLY", "").strip()
TIMEOUT = float(os.environ.get("PROXY_TIMEOUT", "30"))
RETRIES = int(os.environ.get("PROXY_RETRIES", "2"))


def fetch(url, retries=None, verbose=True):
    """依次尝试各代理；返回 (bytes, proxy_name) 或 (None, err)。"""
    if retries is None:
        retries = RETRIES
    enc = urllib.parse.quote(url, safe="")
    last_err = "no proxy configured"
    cand = [(n, t) for (n, t) in PROXIES if (not ONLY or n == ONLY)]
    for attempt in range(retries):
        for name, tmpl in cand:
            tgt = tmpl.format(u=enc, r=url)
            try:
                req = urllib.request.Request(tgt, headers={"User-Agent": UA, "Accept": "*/*"})
                with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                    data = r.read()
                if not data:
                    last_err = f"{name}: empty"
                    continue
                # allorigins /get 返回 JSON 包装
                if name == "allorigins-get":
                    import json
                    try:
                        j = json.loads(data.decode("utf-8", "replace"))
                        c = j.get("contents", "")
                        data = c.encode("utf-8") if isinstance(c, str) else data
                    except Exception as e:
                        last_err = f"{name}: json {e}"
                        continue
                if verbose:
                    print(f"    [ok] {name}  {len(data)} bytes", file=sys.stderr)
                return data, name
            except Exception as e:
                last_err = f"{name}: {e}"
                if verbose:
                    print(f"    [x] {name}: {e}", file=sys.stderr)
        if attempt < retries - 1:
            time.sleep(2 + 2 * attempt)
    return None, last_err


def maybe_gunzip(data, url, outfile):
    """若目标/内容为 gzip 且输出名不带 .gz，则解压。"""
    if outfile.endswith(".gz"):
        return data
    if data[:2] == b"\x1f\x8b":
        try:
            return gzip.decompress(data)
        except Exception:
            return data
    return data


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = {a for a in sys.argv[1:] if a.startswith("--")}
    if not args:
        print(__doc__)
        return 1
    url = args[0]
    if "--list" in flags or len(args) == 1:
        data, info = fetch(url)
        if data is None:
            print("FAIL:", info)
            return 1
        data = maybe_gunzip(data, url, "-")
        sys.stdout.write(data.decode("utf-8", "replace"))
        return 0

    outfile = args[1]
    os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    data, info = fetch(url)
    if data is None:
        print("FAIL:", info)
        return 1
    data = maybe_gunzip(data, url, outfile)
    with open(outfile, "wb") as f:
        f.write(data)
    print(f"OK {outfile}  {len(data)} bytes  via {info}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
