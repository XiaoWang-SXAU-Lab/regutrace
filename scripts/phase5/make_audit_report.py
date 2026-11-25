#!/usr/bin/env python3
"""Generate a self-contained Phase 5 audit HTML (and optional PDF input)."""
import base64
import csv
import html
import json
import os
import statistics

ROOT = os.path.join(os.environ["LOCALAPPDATA"], "rsWGCNA-build", "phase5")
BASE = os.path.join(ROOT, "analysis2")
DAP = os.path.join(ROOT, "dap", "run1")
OUT = os.path.join(os.environ["LOCALAPPDATA"], "rsWGCNA-build", "phase5",
                   "phase5-audit-report.html")
FIG = os.path.join(os.environ["USERPROFILE"], "Bio-info-PRJ", "regutrace",
                   "docs", "development", "figures", "phase4")


def esc(value):
    return html.escape(str(value))


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()


policy = json.load(open(os.path.join(BASE, "phase5-policy-summary.json")))["data"]
audit = read_tsv(os.path.join(BASE, "phase5-coord-audit.tsv"))
targets = read_tsv(os.path.join(BASE, "phase5-policy-targets.tsv"))
null = json.load(open(os.path.join(ROOT, "null", "null-summary.json")))
dbmeta = json.load(open(os.path.join(ROOT, "background", "background-build.json")))
motif = json.load(open(os.path.join(BASE, "phase5-motif-aggregation.json")))
dap_direct = {r["tf_symbol"]: r for r in read_tsv(os.path.join(DAP, "dap-benchmark.tsv")) if r.get("status") == "ok"}
dap_transfer = {r["tf_symbol"]: r for r in read_tsv(os.path.join(DAP, "dap-benchmark-transfer.tsv")) if r.get("status") == "ok"}
dap_emp = {r["tf_symbol"]: r for r in read_tsv(os.path.join(DAP, "dap-benchmark-empirical.tsv")) if r.get("status") == "ok"}

zero = sum(1 for r in audit if r["overlap_threshold"] == "0.8" and int(r["architecture_groups"]) == 0)
one = sum(1 for r in audit if r["overlap_threshold"] == "0.8" and int(r["architecture_groups"]) == 1)


def metric(pol, filt, key):
    return policy[f"{pol}|{filt}"][key]


def pct(value):
    return f"{value * 100:.2f}%"


def bar_rows(items, max_value=None, color="#2563eb", unit=""):
    """items: [(label, value, display)]"""
    maximum = max_value or max(v for _, v, _ in items)
    out = []
    for label, value, display in items:
        width = 0 if maximum == 0 else max(0.0, value / maximum * 100)
        out.append(
            f'<div class="bar-row"><div class="bar-label">{esc(label)}</div>'
            f'<div class="bar-track"><div class="bar-fill" style="width:{width:.1f}%;background:{color}"></div></div>'
            f'<div class="bar-value">{esc(display)}{(" " + esc(unit)) if unit else ""}</div></div>'
        )
    return "\n".join(out)


def grouped_rows(items, colors):
    """items: [(label, [(value, display), ...])]"""
    maximum = max(max(v for v, _ in series) for _, series in items)
    out = []
    for label, series in items:
        bars = "".join(
            f'<div class="bar-fill" style="width:{(v / maximum * 100 if maximum else 0):.1f}%;background:{colors[i]}"></div>'
            for i, (v, _) in enumerate(series)
        )
        labels = " · ".join(f"{esc(d)}" for _, d in series)
        out.append(
            f'<div class="bar-row wide"><div class="bar-label">{esc(label)}</div>'
            f'<div class="bar-track grouped">{bars}</div><div class="bar-value">{labels}</div></div>'
        )
    return "\n".join(out)


# DAP AUPRC grouped bars
dap_names = sorted(set(dap_direct) & set(dap_transfer))
dap_items = []
for name in dap_names:
    dap_items.append((name, [(float(dap_direct[name]["auprc"]), f"{float(dap_direct[name]['auprc']):.3f}"),
                              (float(dap_transfer[name]["auprc"]), f"{float(dap_transfer[name]['auprc']):.3f}")]))

# FDR grouped
fdr_items = []
for pol, label in [("all", "All candidates"), ("best@0.8", "Best-supported"), ("reps@0.8", "Representative groups")]:
    fdr_items.append((label, [(metric(pol, "original", "coverage"), pct(metric(pol, "original", "coverage"))),
                              (metric(pol, "site_fdr", "coverage"), pct(metric(pol, "site_fdr", "coverage"))),
                              (metric(pol, "promoter_fdr", "coverage"), pct(metric(pol, "promoter_fdr", "coverage")))]))

# Fta
reg_all = int(motif["regression_all"]["targets_original"])
reg_best = int(motif["regression_best"]["targets_original"])
reg_reps = int(motif["regression_reps"]["targets_original"])

# DB counts
counts = dbmeta["counts"]
db_rows = "\n".join(
    f'<tr><td><code>{esc(k)}</code></td><td class="num">{v:,}</td></tr>'
    for k, v in counts.items()
)

# Coverage cards
def card(label, value, note):
    return f'<div class="stat"><div class="stat-label">{esc(label)}</div><div class="stat-value">{value}</div><div class="stat-note">{note}</div></div>'

img_overview = b64(os.path.join(FIG, "overview-final.png"))
img_profile = b64(os.path.join(FIG, "tf-profile-final.png"))
FIG5 = os.path.join(os.environ["USERPROFILE"], "Bio-info-PRJ", "regutrace", "docs", "development", "figures", "phase5")
img_policy = b64(os.path.join(FIG5, "policy-ui-best.png"))

html_doc = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>ReguTrace Phase 5 审计报告</title>
<style>
:root {{
  --ink:#172033; --muted:#5f6b7a; --line:#dce3ec; --bg:#f5f7fb; --card:#ffffff;
  --blue:#2563eb; --orange:#d97706; --teal:#0f766e; --red:#b91c1c; --green:#15803d;
}}
* {{ box-sizing: border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink);
  font-family:"Microsoft YaHei","PingFang SC","Noto Sans CJK SC","Segoe UI",Arial,sans-serif;
  font-size:14px; line-height:1.65; }}
.page {{ max-width:1120px; margin:0 auto; padding:34px 38px 60px; }}
.hero {{ background:linear-gradient(135deg,#0f2742,#164e63); color:#fff; border-radius:14px;
  padding:30px 34px; margin-bottom:22px; }}
.hero h1 {{ margin:0 0 6px; font-size:30px; letter-spacing:.02em; }}
.hero .sub {{ color:#c7d7e8; font-size:15px; }}
.hero .meta {{ margin-top:14px; color:#dbeafe; font-size:12px; }}
h2 {{ font-size:21px; margin:30px 0 10px; border-left:5px solid var(--blue); padding-left:11px; }}
h3 {{ font-size:16px; margin:20px 0 8px; }}
p {{ margin:7px 0; }}
.small {{ font-size:12px; color:var(--muted); }}
.muted {{ color:var(--muted); }}
.grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; }}
.grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:10px; padding:16px 18px; margin:10px 0; }}
.stat {{ background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }}
.stat-label {{ font-size:12px; color:var(--muted); }}
.stat-value {{ font-size:25px; font-weight:700; margin:2px 0; }}
.stat-note {{ font-size:11px; color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; margin:8px 0 4px; font-size:12.5px; }}
th,td {{ border-bottom:1px solid var(--line); padding:7px 8px; vertical-align:top; text-align:left; }}
th {{ background:#eef2f7; font-weight:600; }}
td.num {{ text-align:right; font-variant-numeric:tabular-nums; white-space:nowrap; }}
code {{ background:#eef2f7; padding:1px 4px; border-radius:4px; font-family:ui-monospace,Consolas,monospace; font-size:11.5px; }}
.badge {{ display:inline-block; padding:1px 7px; border-radius:999px; font-size:11px; font-weight:600; }}
.ok {{ background:#dcfce7; color:#166534; }}
.warn {{ background:#fef3c7; color:#92400e; }}
.bad {{ background:#fee2e2; color:#991b1b; }}
.info {{ background:#dbeafe; color:#1e40af; }}
.bar-row {{ display:grid; grid-template-columns:150px 1fr 110px; gap:10px; align-items:center; margin:7px 0; }}
.bar-row.wide {{ grid-template-columns:150px 1fr 230px; }}
.bar-label {{ font-size:12px; color:#334155; }}
.bar-track {{ height:16px; background:#e8edf4; border-radius:4px; overflow:hidden; display:flex; }}
.bar-track.grouped {{ gap:2px; }}
.bar-fill {{ height:100%; border-radius:3px; }}
.bar-value {{ font-size:11px; color:var(--muted); font-variant-numeric:tabular-nums; }}
.figure {{ background:#fff; border:1px solid var(--line); border-radius:10px; padding:14px 16px; margin:12px 0; }}
.figure .cap {{ font-size:11.5px; color:var(--muted); margin-top:6px; }}
.flow {{ display:flex; gap:8px; align-items:stretch; flex-wrap:wrap; }}
.flow .node {{ flex:1; min-width:140px; background:#eef6ff; border:1px solid #bfdbfe; border-radius:8px; padding:10px; }}
.flow .arrow {{ align-self:center; color:#64748b; font-size:20px; }}
.flow .node b {{ display:block; color:#1e40af; }}
.callout {{ border-left:4px solid var(--orange); background:#fffbeb; padding:10px 14px; border-radius:0 8px 8px 0; margin:10px 0; }}
.danger {{ border-left-color:var(--red); background:#fef2f2; }}
.good {{ border-left-color:var(--green); background:#f0fdf4; }}
ul {{ margin:6px 0 6px 20px; padding:0; }}
li {{ margin:3px 0; }}
.two-col {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
img.shot {{ width:100%; border:1px solid var(--line); border-radius:8px; }}
.footer {{ margin-top:34px; padding-top:12px; border-top:1px solid var(--line); color:var(--muted); font-size:11px; }}
@media print {{
  body {{ background:#fff; font-size:11px; }}
  .page {{ max-width:none; padding:0; }}
  .hero {{ break-inside:avoid; }}
  h2 {{ break-after:avoid; }}
  .card,.figure,.stat {{ break-inside:avoid; }}
  @page {{ margin:14mm; }}
}}
</style>
</head>
<body><div class="page">
<div class="hero">
  <h1>ReguTrace Phase 5 审计报告</h1>
  <div class="sub">科学声明、代码、数据与可复用 regulatory background 的逐项复核</div>
  <div class="meta">审计日期：2026-11-25 · 代码：regutrace 67379cf + Final Closure · dashboard f23f766 + Final Closure<br>
  性质：开发代理自审（非独立第三方）；供 ChatGPT 二次 review 与后续 Phase 6 决策</div>
</div>

<div class="callout good"><b>审计结论：PASS WITH RESERVATIONS（有条件通过）。</b>
Phase 5 的主要数值、target policy 对比、可复用 SQLite background 和 round-trip 工程链路均可复现；
但 BH q 值是条件性上界而非无条件 FDR 控制、DAP benchmark 没有真实负例、坐标覆盖仍有 52/1,143 个 TF 缺口；
v3 已精确保存 all/best/reps 并通过 per-TF 集合核验，完整项目导入与重启/reopen 已实测通过。当前结果足以支持“独立候选 regulatory background”，
不足以支持“校准后的体内调控网络”或“全基因组精度”声明。</div>

<h2>1. 审计范围与方法</h2>
<div class="card">
<p><b>审阅对象：</b>Phase 5 报告及 supporting 文档、<code>scripts/phase5/</code> 的分析/导出/验证脚本、
HF1 background SQLite、Phase 4 dashboard 集成、基线提交 <code>67379cf</code> 与 <code>f23f766</code>，
以及 Final Closure 改动（审计更正、统计测试、UI motif policy）。</p>
<p><b>审计方法：</b>①静态阅读代码与文档；②不信任摘要，直接从 <code>phase5-policy-targets.tsv</code> 重算覆盖和分位数；
③从 site table 重算 FDR 目标集合的抽样核对；④直接查询 SQLite 的 schema、行数、hash 与索引；
⑤实测 <code>validate-db</code> 和 <code>regutrace prior</code>；⑥对照 DAP-seq 明细重算中位数；
⑦检查 macOS 路径测试修复的可移植性。</p>
<p><b>未执行：</b>macOS 原生构建/UI 验收、FIMO 实跑、desktop WebView 对 policy 选择器的截图验收、
全基因组 dependence-aware null。</p>
</div>

<h2>2. 数据流与证据边界</h2>
<div class="figure">
<div class="flow">
  <div class="node"><b>Genome + GFF</b>37,415 promoters<br>1,000 bp upstream</div><div class="arrow">→</div>
  <div class="node"><b>TF catalogue</b>1,981 rows<br>1,907 identified</div><div class="arrow">→</div>
  <div class="node"><b>TF–motif map</b>32,792 rows<br>strict/calibrated/exploratory</div><div class="arrow">→</div>
  <div class="node"><b>Site scan</b>8,191,211 unique sites<br>677 PWM groups</div><div class="arrow">→</div>
  <div class="node"><b>Background</b>1.41 GB SQLite<br>4.64M motif–target / 7.39M pairs<br>44,295 policy rows</div>
</div>
<div class="cap">图 1：Phase 5 复用的数据流。所有 PWM 派生 target 关系仍为 evidence level D；map 的 strict/calibrated/exploratory 状态与 target 证据等级分开。</div>
</div>

<div class="callout"><b>关键边界：</b>本次审计把“TF 鉴定”“motif 映射”“可测试性”“候选 target”视为四个独立层；
<code>58.94%</code>、<code>51.97%</code>、<code>57.89%</code> 都是 actual target coverage，不是 precision、recall 或实验验证率。</div>

<div class="card">
<h3>2.1 四个计数层（不得混用）</h3>
<table><thead><tr><th>层</th><th>数量</th><th>含义</th></tr></thead><tbody>
<tr><td>TF classification</td><td class="num">1,907 identified (+74 ambiguous)</td><td>catalogue 级身份；不等于有 motif</td></tr>
<tr><td>Motif mapping</td><td class="num">1,143 TFs；679 mapped motifs</td><td>至少一个受支持 motif 分配</td></tr>
<tr><td>Coordinate availability</td><td class="num">1,089 joined；2 key mismatch；52 unavailable</td><td>audit 的 1,091 个 one-group TF 不等于 DB join 数；22 个 coord-only TF 不在 baseline</td></tr>
<tr><td>Actual target coverage</td><td class="num">1,124 / 1,907</td><td>p≤1e-4 下至少一个 promoter candidate</td></tr>
</tbody></table>
<p class="small">坐标状态应拆成：<b>1,089 successfully joined + 2 key mismatches + 52 unavailable</b>。1,091 是 audit 中有 1 个 DBD group 的 base TF 数，不等同于 DB join 数；22 个 coord-only TF 不在 baseline。52 个无坐标 base TF 仍计入 1,143 mapping/coverage 分母，但在 best/reps 中使用 opaque fallback。</p>
<p class="small"><code>domain_mode=single_domain_fallback</code> 的 assignment 是完整 DBD 架构提取失败后的 alignment-based fallback；UI 已标注，且不会按 full architecture 参与排序或升级。</p>
</div>

<h2>3. 声明逐项核验</h2>
<div class="card">
<table>
<thead><tr><th>声明</th><th>来源</th><th>独立核验</th><th>判定</th></tr></thead>
<tbody>
<tr><td>HF1 分母 1,907；CDS 1,908 只差 <code>Fta_6G0009940</code></td><td>phase5-report §1</td><td>比较 CDS/PEP catalogue 的 identified 集合；唯一差异确认为该基因（CDS AP2；PEP ambiguous AP2;ERF）</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>All policy：1,124 TFs / 58.94% / 25,481,360 pairs</td><td>policy summary</td><td>从 <code>phase5-policy-targets.tsv</code> 对每个 TF 的 target 数独立求和；完全一致</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>Best@0.8：991 / 51.97% / 7,394,435</td><td>policy summary</td><td>独立重算；中位 4,091，Q1/Q3 2,290/8,385；一致</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>Reps@0.8：1,104 / 57.89% / 19,848,826</td><td>policy summary</td><td>独立重算；中位 14,657；一致</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>v3 policy exactness：all/best/reps 的 TF 集合与 per-TF target 集合</td><td>policy verification</td><td>从 phase5-policy-targets.tsv 读取选择，按 TF 与 DB join 逐集合比较；0 mismatch；pair-set hash 已记录</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>FDR site/promoter：All 22.50/21.87；Best 4.51/4.25；Reps 17.09/15.05</td><td>policy summary</td><td>独立重算 target 集合与 summary 一致；q 实现见第 6 节；是条件性上界，不是无条件 FDR 控制</td><td><span class="badge warn">PASS WITH CAVEAT</span></td></tr>
<tr><td>Coordinate：每个有坐标 TF 只有 1 个 DBD region group</td><td>coord audit</td><td><b>1,089 joined + 2 key mismatch + 52 unavailable</b>；1,091 是 audit one-group 数；22 个 coord-only TF 不在 baseline</td><td><span class="badge warn">PASS WITH CAVEAT</span></td></tr>
<tr><td><code>Fta_1G0005130</code>：19 motifs / 7 specificity groups / 8,709→2,290→4,112</td><td>motif aggregation</td><td>从 policy target 行与 coord audit 独立核对；pairwise Jaccard 中位 0.0018</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>DAP direct 中位 AUPRC 0.330；transfer 0.261；empirical 0.346</td><td>DAP report</td><td>从三张 per-TF TSV 重算中位数；一致；7 TFs，无真实负例</td><td><span class="badge warn">PASS WITH CAVEAT</span></td></tr>
<tr><td>Background DB 1.41 GB (v3) / validate-db / prior query</td><td>roundtrip</td><td>实测 hash、行数、validate-db、prior、integrity_check</td><td><span class="badge ok">PASS</span></td></tr>
<tr><td>项目导入 / 重启 / upstream 查询</td><td>roundtrip</td><td>v3 DB 实测：1.4 GB import 74.4 s、restart/reopen、753 upstream rows；输入 hash = 7a33…</td><td><span class="badge ok">PASS</span></td></tr>
</tbody>
</table>
</div>

<h2>4. HF1 三种 policy 的实际结果</h2>
<div class="grid">
{card("All candidates", "58.94%", "1,124 / 1,907 TFs · 25.48M pairs")}
{card("Best-supported", "51.97%", "991 / 1,907 TFs · 7.39M pairs")}
{card("Representative groups", "57.89%", "1,104 / 1,907 TFs · 19.85M pairs")}
{card("Median target reduction", "81.16%", "21,715 → 4,091; 18.84% remains")}
</div>
<div class="figure">
<h3>4.1 覆盖与 target 数（原始 p≤1e-4）</h3>
{bar_rows([
  ("All candidates", metric("all", "original", "coverage"), pct(metric("all", "original", "coverage"))),
  ("Best-supported", metric("best@0.8", "original", "coverage"), pct(metric("best@0.8", "original", "coverage"))),
  ("Representative groups", metric("reps@0.8", "original", "coverage"), pct(metric("reps@0.8", "original", "coverage"))),
], 0.65, "#2563eb")}
<div class="cap">图 2：actual target coverage（分母 1,907 identified TFs）。Best 将 pair 总量压缩 71.0%，仍保留 51.97% 覆盖；这不是“减少 target 越多越好”的优化，而是可解释的假设选择。</div>
</div>

<div class="card">
<table>
<thead><tr><th>Policy</th><th>≥1 target</th><th>Coverage</th><th>Median motifs/TF</th><th>Median targets/TF</th><th>Q1 / Q3</th><th>Pairs</th></tr></thead>
<tbody>
<tr><td>All candidates</td><td class="num">1,124</td><td class="num">58.94%</td><td class="num">17</td><td class="num">21,715</td><td class="num">13,891 / 34,348</td><td class="num">25,481,360</td></tr>
<tr><td>Best@0.8</td><td class="num">991</td><td class="num">51.97%</td><td class="num">1</td><td class="num">4,091</td><td class="num">2,290 / 8,385</td><td class="num">7,394,435</td></tr>
<tr><td>Reps@0.8</td><td class="num">1,104</td><td class="num">57.89%</td><td class="num">6</td><td class="num">14,657</td><td class="num">7,644 / 28,711</td><td class="num">19,848,826</td></tr>
</tbody>
</table>
</div>
<div class="callout"><b>审计判断：</b>policy 差异主要来自同一 DBD 区的 donor/PWM 冗余，不是多个独立 DBD。Best 不是“真实供体”的生物学断言；它是可复现的序列/比对排序，并保留全部原始候选供审计。</div>

<h2>5. 坐标审计与回归案例</h2>
<div class="figure">
<h3>5.1 Query DBD group 数（1,143 base TFs，overlap 0.8）</h3>
{bar_rows([
  ("1 个 group", one, f"{one:,}"),
  ("无坐标行", zero, f"{zero:,}"),
], 1200, "#0f766e")}
<div class="cap">图 3：1,091 个有坐标 TF 均为 1 个 query DBD region；52 个 base TF 没有坐标行。坐标阈值 0.5/0.7/0.8/0.9 在 HF1 上不改变分组。52 个缺口会使 best/reps 对这些 TF 退化为“opaque single group”，见 P2-1。</div>
</div>

<div class="grid2">
<div class="figure">
<h3>5.2 Fta_1G0005130：19 个 motif 的 target union</h3>
{bar_rows([
  ("All (19 motifs)", reg_all, f"{reg_all:,}"),
  ("Best (1 motif)", reg_best, f"{reg_best:,}"),
  ("Reps (7 groups)", reg_reps, f"{reg_reps:,}"),
], reg_all, "#d97706")}
<div class="cap">同一 bHLH DBD（152–198），19 个不同 PWM，7 个 specificity groups。target 集合几乎不重叠：median pairwise Jaccard {motif['regression_pairwise_jaccard_median']:.4f}。</div>
</div>
<div class="card">
<h3>5.3 审计发现：坐标覆盖不是 30</h3>
<p>独立核对 <code>phase5-coord-audit.tsv</code> 得到 <b>52</b> 个 base TFs 的 <code>architecture_groups=0</code>，而非报告早期写的 30。
另有 <b>22</b> 个 coordinate-rerun TFs 不在 1,143-row baseline 内。坐标状态为 <b>1,089 successfully joined + 2 key mismatches + 52 unavailable</b>；1,091 是 audit 的 one-group TF 数。2 个 mismatch 是 <code>Fta_4G0031600</code>、<code>Fta_5G0037620</code>（坐标行与 base 的 motif/reference key 不一致）。</p>
<p><b>影响：</b>52 个 unavailable TF 的 best/reps 不能用坐标分组，使用 documented opaque single-group fallback（best/reps 各选一个最佳 exploratory motif）；它们仍被计入 1,143 和 target 统计。这是 <span class="badge warn">P2</span>，已在文档中纠正。</p>
</div>
</div>

<h2>6. Multiple testing 与 promoter-level 统计</h2>
<div class="figure">
<h3>6.1 原始 p 与 BH q≤0.05 后的 coverage</h3>
{grouped_rows(fdr_items, ["#2563eb", "#d97706", "#0f766e"])}
<div class="cap">图 4：蓝=原始 p≤1e-4；橙=site BH q≤0.05；绿=promoter BH q≤0.05。完整测试 universe 为 37,415 promoters × 双链 × 每 motif 1,000 bp windows；10 bp motif 约 74.16M windows。</div>
</div>

<div class="card">
<table>
<thead><tr><th>Policy</th><th>Original coverage</th><th>Site BH coverage</th><th>Promoter BH coverage</th><th>Original pairs</th><th>Site BH q≤0.05 pairs</th><th>Promoter BH q≤0.05 pairs</th></tr></thead>
<tbody>
<tr><td>All</td><td class="num">58.94%</td><td class="num">22.50%</td><td class="num">21.87%</td><td class="num">25,481,360</td><td class="num">2,842,315</td><td class="num">2,174,252</td></tr>
<tr><td>Best</td><td class="num">51.97%</td><td class="num">4.51%</td><td class="num">4.25%</td><td class="num">7,394,435</td><td class="num">393,946</td><td class="num">264,936</td></tr>
<tr><td>Reps</td><td class="num">57.89%</td><td class="num">17.09%</td><td class="num">15.05%</td><td class="num">19,848,826</td><td class="num">1,293,478</td><td class="num">945,875</td></tr>
</tbody>
</table>
</div>

<div class="grid2">
<div class="figure">
<h3>6.2 Dinucleotide shuffle null</h3>
{bar_rows([
  ("Observed promoters", null["observed_hit_rate"], f"{null['observed_hit_rate']*100:.1f}%"),
  ("Shuffled promoters", null["null_hit_rate"], f"{null['null_hit_rate']*100:.1f}%"),
], 0.7, "#7c3aed")}
<div class="cap">200 promoters，3 个 motif，200 kb。观察 117/200，shuffle 123/200；观察组并不高于 composition-matched 序列。该 null 是示例性，不是全基因组零模型。</div>
</div>
<div class="card">
<h3>6.3 审计判断</h3>
<ul>
<li><span class="badge warn">P1（科学）</span> site/promoter q 是条件性上界：只对 hit 排序；在 hit 表完整、window 数精确且无缺失检验时可作为 full BH 的上界，但不建模窗口/基因依赖、GC/长度或 motif 家族依赖，也不提供无条件 FDR 控制。</li>
<li><span class="badge ok">已披露</span> 默认 p=1e-4 未变；BH q≤0.05 仅作为可选解释层，不能替代原始候选背景。</li>
<li><span class="badge warn">缺口</span> 没有 block/bootstrap 或 dependence-aware null；FIMO 未安装/未执行。</li>
</ul>
</div>
</div>

<h2>7. 独立 DAP-seq benchmark</h2>
<div class="figure">
<h3>7.1 Direct vs held-out DBD transfer（AUPRC，6 个共同 TF）</h3>
{grouped_rows(dap_items, ["#2563eb", "#d97706"])}
<div class="cap">图 5：蓝=direct JASPAR motif；橙=held-out DBD transfer（benchmark TFs 从 donor map 移除）。图为 6 个共同 TF；RAV1 只有 transfer，GBF3 只有 direct。Direct 中位 0.330；transfer 0.261；TGA1 transfer 反而更好（0.449 vs 0.255）。</div>
</div>

<div class="grid">
{card("Direct median AUPRC", "0.330", "7 TFs, uniform")}
{card("Transfer median AUPRC", "0.261", "7 TFs, held-out")}
{card("Empirical AUPRC", "0.346", "higher recall, lower PPV")}
{card("Common-TF direct", "0.314", "vs transfer 0.279")}
</div>

<div class="callout danger"><b>P1（科学限制）：</b>DAP-seq 的 positive 是 promoter 与 narrowPeak 重叠；未观测基因是 unknown，不是 negative。AUPRC/PPV_assay 使用 unknown-as-non-positive；偏差方向未建立（未观测可能是假阴性，assay peak 也可能是间接结合），因此只能作 assay-defined 比较，不能称为生物学 precision。只覆盖 7 个 TF、一个 assay series，且未跑 Arabidopsis strict/calibrated arm。</div>

<h2>8. 可复用 background 与 round-trip</h2>
<div class="grid2">
<div class="card">
<h3>8.1 SQLite schema</h3>
<table><tbody>{db_rows}</tbody></table>
<p class="small">schema v3；<code>hits=0</code> 是 site-only 设计，<code>pairs</code> 是 best policy，不是 25M union。<code>tf_motif_policy</code> 保存 all/best/reps 的 44,295 条精确 TF–motif 选择；<code>tf_motif_coordinate</code> 保存 28,465 条 coordinate rerun；<code>tf_motif_policy_config</code> 记录 policy version、参数和源文件 hash。目标集由 join <code>motif_target</code> 重建，不物化三份 TF–target 表。精确核验：all 1,124 / 25,481,360；best 991 / 7,394,435；reps 1,104 / 19,848,826，per-TF 目标集合 0 mismatch。</p>
</div>
<div class="card">
<h3>8.2 实测结果</h3>
<table><tbody>
<tr><th>大小</th><td class="num">{dbmeta['database_bytes']:,} bytes</td></tr>
<tr><th>SHA-256</th><td><code>{esc(dbmeta.get('database_sha256',''))}</code></td></tr>
<tr><th>validate-db</th><td><span class="badge ok">PASS</span> (v3；integrity_check ok)</td></tr>
<tr><th>独立 prior query</th><td>3 genes → 321 rows, 1.03 s</td></tr>
<tr><th>项目 Data Import (v3)</th><td><span class="badge ok">PASS</span> 74.4 s</td></tr>
<tr><th>重启后 project upstream</th><td><span class="badge ok">200</span> 753 rows, 2.1 s</td></tr>
<tr><th>第二 independent project</th><td><span class="badge ok">200</span> 753 rows, 2.1 s</td></tr>
<tr><th>SHA-256 一致性</th><td>项目输入文件 = <code>7a33d3fc…</code></td></tr>
</tbody></table>
</div>
</div>

<div class="callout danger"><b>审计发现（已修复）：</b>原 <code>motif_registry</code> 只有 805 行，来自 JASPAR 2024；<code>tf_motif</code> 的 679 个 mapped motif 中有 111 个（2026 ID）没有 registry 元数据。已用实际 927-row scan library 重建 registry，现 0 个 mapped motif 缺失；最终 v3 hash 为 <code>7a33d3fc…</code>。这是本次审计最重要的可修复数据缺口。</div>

<h2>9. 测试与可复现性</h2>
<table>
<thead><tr><th>检查</th><th>结果</th><th>备注</th></tr></thead>
<tbody>
<tr><td>Rust <code>cargo test --workspace</code></td><td><span class="badge ok">124 passed</span></td><td>46 CLI unit + 7 CLI + 2 Phase 3 + 69 core</td></tr>
<tr><td>Python Phase 3</td><td><span class="badge ok">6 passed</span></td><td>synthetic validation</td></tr>
<tr><td>Phase 5 unit</td><td><span class="badge ok">5 passed</span></td><td>grouping/overlap/dedup/policy</td></tr>
<tr><td>Phase 5 statistics</td><td><span class="badge ok">7 passed</span></td><td>complete universe/ties/monotonic/promoter</td></tr>
<tr><td>Go regutrace</td><td><span class="badge ok">56 passed</span></td><td>Windows；含 10 个 motif-policy 测试</td></tr>
<tr><td>Go api</td><td><span class="badge ok">13 passed</span></td><td>Windows</td></tr>
<tr><td>Frontend build</td><td><span class="badge ok">pass</span></td><td>vue-tsc；vite build 见 CLI</td></tr>
<tr><td>macOS path fix</td><td><span class="badge warn">compiled on Windows</span></td><td>f23f766；未在 Mac 运行</td></tr>
<tr><td>v3 DB import</td><td><span class="badge ok">pass</span></td><td>1.4 GB import；重启/reopen 后查询成功</td></tr>
<tr><td>synthetic UI policy acceptance</td><td><span class="badge ok">all 6 / best 4 / reps 5</span></td><td>实际 release dashboard 切换 selector；fixed strict/calibrated 保留，ambiguous-only G6 在 best/reps 排除</td></tr>
</tbody>
</table>

<div class="two-col">
<div class="figure">
<h3>UI 证据（Phase 4，非 Phase 5 新 UI）</h3>
<img class="shot" src="data:image/png;base64,{img_overview}" alt="ReguTrace overview">
<div class="cap">概览：coverage layers、mapping status、background 和 site-only 标签。</div>
</div>
<div class="figure">
<h3>TF profile（Phase 4）</h3>
<img class="shot" src="data:image/png;base64,{img_profile}" alt="ReguTrace TF profile">
<div class="cap">真实 PWM logo、testability、shared-PWM 警告和分页 target 列表。此图为加入 policy 选择器前的 Phase 4 界面。</div>
</div>
</div>
<div class="figure">
<h3>Final Closure：生产 policy 选择器（synthetic runtime acceptance）</h3>
<img class="shot" src="data:image/png;base64,{img_policy}" alt="ReguTrace policy selector">
<div class="cap">Best-supported 选中后：target 4（strict/calibrated G1/G2 + exploratory G3/G5），policy 信息显示 4 motif / 2 fixed / 2 exploratory；all=6、reps=5。截图来自实际运行的 release dashboard；原始 site p≤1e-4 与全量 evidence 导出仍不按 policy 过滤。</div>
</div>

<h2>10. 审计发现与严重度</h2>
<div class="card">
<table>
<thead><tr><th>ID</th><th>严重度</th><th>发现</th><th>影响/建议</th></tr></thead>
<tbody>
<tr><td>P1-1</td><td><span class="badge bad">P1</span></td><td>DAP benchmark 无真实负例、仅 promoter overlap、7 TFs，且未跑 strict/calibrated arm</td><td>只能作 comparative diagnostic；论文中不能写“precision”；Phase 6 前扩展 ChIP/DAP 与距离规则</td></tr>
<tr><td>P1-2</td><td><span class="badge bad">P1</span></td><td>site/promoter BH 是条件性上界；未建模 dependence/GC/length</td><td>保留为可选层；已有 7 项统计回归测试；block bootstrap 与 GC/length-matched null 仍未做</td></tr>
<tr><td>P2-1</td><td><span class="badge warn">P2</span></td><td>坐标状态 = 1,089 joined + 2 key mismatch + 52 unavailable；22 个 coord-only TFs 不在 baseline</td><td>best/reps 对 52 个 unavailable TF 退化为 opaque single group；需补坐标或明确排除</td></tr>
<tr><td>P2-2</td><td><span class="badge warn">P2</span></td><td>ranking 未使用 full-length similarity、reciprocity、DNA-contact conservation（字段不可得/未实现）</td><td>当前 best 仅是 architecture/coverage/identity/bitscore 排序；不能声称使用了全部 requested ranking factors</td></tr>
<tr><td>P2-3</td><td><span class="badge ok">已修复</span></td><td>v3 DB 完整 dashboard 项目导入、重启/reopen 已执行</td><td>输入 SHA-256 = 7a33…；查询 753 rows</td></tr>
<tr><td>P2-4</td><td><span class="badge warn">P2</span></td><td>ReguTrace job registry 仍为内存态；job ID 重启不恢复</td><td>导出 DB 已可复用；如需 job 重开，需持久化+懒加载</td></tr>
<tr><td>P2-5</td><td><span class="badge warn">P2</span></td><td>旧 <code>motif_registry</code> 只有 805 行，111 个 mapped 2026 motif 无元数据</td><td><b>已修复</b>：用 927-row job library 重建；现 0 missing</td></tr>
<tr><td>P3-1</td><td><span class="badge ok">已修复</span></td><td>layers/policy 原先需另跑脚本</td><td><code>build_background_v2.py</code> 串联四步并记录最终 hash/counts</td></tr>
<tr><td>P3-2</td><td><span class="badge info">P3</span></td><td>macOS 测试修复只在 Windows 编译验证</td><td>Mac 端跑 <code>go test ./internal/regutrace</code></td></tr>
<tr><td>P3-3</td><td><span class="badge info">P3</span></td><td>JASPAR 2026 默认化未决；FIMO 未安装</td><td>保持 2024 默认；Phase 6 或专门 benchmark 再决定</td></tr>
</tbody>
</table>
</div>

<h2>11. 给 ChatGPT 的复核问题</h2>
<div class="card">
<ol>
<li>保守 BH 的公式（hit ranks + complete window universe）是否可作为 Phase 5 的可接受上界？是否应改为 block/bootstrap 或直接报告未校正候选？</li>
<li>DAP-seq 的 positive=promoter overlap、unknown=unobserved 是否足以支持“direct motif 优于 transfer”的相对结论？还缺哪些对照？</li>
<li>52 个无坐标 TF 的 fallback 处理是否应排除在 best/reps 覆盖之外，还是保留但仍算 coverage？</li>
<li><code>pairs</code> 仅导出 best policy（7.39M），而 <code>motif_target</code> 保留 4.64M motif 关系，是否满足 rsWGCNA 的“保留原始证据、下游自行过滤”要求？</li>
<li>v3 的 <code>tf_motif_policy</code> + <code>tf_motif_coordinate</code> 是否足以作为 Phase 6 的唯一 policy 接口，还是还需要在导出层提供按 policy 物化的 TSV？</li>
<li>Phase 6 的 rsWGCNA integration 是否应先解决 dependence-aware statistics，还是先做 expression-only GRN 的证据分层？</li>
</ol>
</div>

<h2>12. 复现命令与文件位置</h2>
<div class="card">
<pre class="small">regutrace validate-db --db analysis/regutrace-hf1-phase5-background-2026-11-25/hf1-background.db
regutrace prior --db ... --targets targets-example.txt --out upstream.tsv

python scripts/phase5/aggregate_audit.py --coord-map ... --base-map ... --sites ... \
  --testability ... --pwm-hash ... --specificity-groups ... --out-dir ...

python scripts/phase5/build_background_v2.py ... --policy best@0.8 \
  --promoter-summary ... --scan-manifest ... --motifs <job motifs.jaspar> \
  --policy-targets <phase5-policy-targets.tsv> --policy-version phase5-v1
# wrapper runs build_background.py -> add_background_layers.py -> rebuild_motif_registry.py -> add_policy_layers.py

python scripts/phase5/verify_policy_layers.py --db hf1-background.db \
  --policy-targets phase5-policy-targets.tsv \
  --summary phase5-policy-summary.json

cargo test --workspace
go test ./internal/regutrace ./internal/api
pnpm build</pre>
<p class="small">阶段基线：regutrace <code>67379cf</code>，dashboard <code>f23f766</code>；Final Closure 的审计更正、统计测试与 UI motif policy 已纳入本轮 checkpoint。Background SHA-256 <code>{esc(dbmeta.get('database_sha256',''))}</code>。</p>
</div>

<div class="card"><h3>Closure 后仍需完成</h3>
<ul>
<li><b>Manuscript-level：</b>更大规模独立 DAP/ChIP benchmark、FIMO 实跑、dependence-aware promoter 统计、strict/calibrated Arabidopsis arm、reps target-set 外部验证。</li>
<li><b>macOS native acceptance：</b>macOS `go test`、`/var` 路径回归、WebView policy selector/fallback/import-reopen 验收；本轮未重建安装器。</li>
<li><b>Schema freeze：</b>v3（policy version `phase5-v1`）；后续修改需新 policy/schema version，不得静默改变已有 exact 表。</li>
</ul></div>

<div class="footer">
本报告是开发代理自审，不是独立第三方审计，也不是实验验证。所有 PWM motif–target 关系为 evidence level D；
未观测基因不是实验阴性；未做 FIMO 实跑、macOS 原生 UI 验收或 Phase 6 网络分析。
</div>
</div></body></html>
"""

open(OUT, "w", encoding="utf-8").write(html_doc)
print(OUT, len(html_doc))
