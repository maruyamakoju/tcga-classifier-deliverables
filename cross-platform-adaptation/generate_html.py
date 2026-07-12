#!/usr/bin/env python
"""Generate a self-contained, JS-free HTML summary (inline SVG) from the
benchmark CSVs. ASCII-only output so it renders identically from file://."""
import os
import csv
import html

HERE = os.path.dirname(os.path.abspath(__file__))
def read_csv(name):
    with open(os.path.join(HERE, name), newline="") as f:
        return list(csv.DictReader(f))

bench = read_csv("adaptation_benchmark.csv")
imb = read_csv("adaptation_imbalance.csv")

INK = "#1f2933"; MUT = "#8a94a6"; GRID = "#e3e8ef"
DEPLOY = "#9aa5b1"      # muted grey  = deployed / no adaptation
ADAPT = "#2f6df6"       # strong blue = cohort standardization
GOOD = "#12b886"; BAD = "#e8590c"

def bar_group(title, subtitle, cats, series, ymax=1.0, fmt="{:.3f}", unit=""):
    """series: list of (label, color, [values aligned to cats])."""
    W, H = 560, 300
    padL, padR, padT, padB = 48, 16, 46, 58
    plotW, plotH = W - padL - padR, H - padT - padB
    n = len(cats); gw = plotW / n
    nb = len(series); bw = gw * 0.62 / nb
    def y(v): return padT + plotH * (1 - v / ymax)
    svg = [f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="{html.escape(title)}">']
    svg.append(f'<text x="{padL}" y="20" class="ttl">{html.escape(title)}</text>')
    svg.append(f'<text x="{padL}" y="36" class="sub">{html.escape(subtitle)}</text>')
    for gv in [0, .25, .5, .75, 1.0]:
        yy = y(gv * ymax)
        svg.append(f'<line x1="{padL}" y1="{yy:.1f}" x2="{W-padR}" y2="{yy:.1f}" class="grid"/>')
        svg.append(f'<text x="{padL-6}" y="{yy+3:.1f}" class="axl" text-anchor="end">{gv*ymax:.2f}</text>')
    for ci, cat in enumerate(cats):
        gx = padL + ci * gw
        for si, (lab, col, vals) in enumerate(series):
            v = vals[ci]
            if v is None: continue
            bx = gx + gw * 0.19 + si * bw
            by = y(v); bh = padT + plotH - by
            svg.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw*0.86:.1f}" height="{bh:.1f}" rx="2" fill="{col}"/>')
            svg.append(f'<text x="{bx+bw*0.43:.1f}" y="{by-4:.1f}" class="val" text-anchor="middle">{fmt.format(v)}{unit}</text>')
        svg.append(f'<text x="{gx+gw/2:.1f}" y="{H-padB+16:.1f}" class="cat" text-anchor="middle">{html.escape(cat)}</text>')
    # legend
    lx = padL
    for lab, col, _ in series:
        svg.append(f'<rect x="{lx}" y="{H-18}" width="11" height="11" rx="2" fill="{col}"/>')
        svg.append(f'<text x="{lx+16}" y="{H-8}" class="leg">{html.escape(lab)}</text>')
        lx += 30 + len(lab) * 7
    svg.append('</svg>')
    return "\n".join(svg)

# ---- data for charts ----
def getrow(cohort_sub, mode_sub):
    for r in bench:
        if cohort_sub in r["cohort"] and mode_sub in r["mode"]:
            return r
    return None

# Chart 1: Toil acc & specificity, deployed vs adapted
toil_dep = getrow("Toil", "deployed"); toil_da = getrow("Toil", "z-score")
c1 = bar_group(
    "Historical TCGA-Toil/RSEM cohort: threshold behavior after adaptation",
    "n=200, 100 tumor / 100 normal. Observed AUC 0.992 -> 0.994.",
    ["Accuracy @0.5", "Balanced acc", "1 - FPR (specificity)"],
    [("Deployed (no adaptation)", DEPLOY,
      [float(toil_dep["acc_at_0p5"]), float(toil_dep["balanced_acc"]), 1-float(toil_dep["FPR_at_0p5"])]),
     ("Cohort standardization", ADAPT,
      [float(toil_da["acc_at_0p5"]), float(toil_da["balanced_acc"]), 1-float(toil_da["FPR_at_0p5"])])])

# Chart 2: imbalance robustness (balanced acc vs tumor fraction)
fracs = [r["tumor_fraction"] for r in imb]
c2 = bar_group(
    "Threshold recovery depends on cohort composition",
    "Toil resampled to each tumor fraction; balanced accuracy at 0.5 (AUC stays ~0.995 throughout).",
    [f"{float(x)*100:.0f}% tumor" for x in fracs],
    [("Deployed", DEPLOY, [float(r["balanced_acc_deployed"]) for r in imb]),
     ("Cohort standardization", ADAPT, [float(r["balanced_acc_adapted"]) for r in imb])])

# GTEx stat
gtex_dep = getrow("GTEx", "deployed"); gtex_da = getrow("GTEx", "z-score")
gtex_fpr_dep = float(gtex_dep["FPR_at_0p5"]); gtex_fpr_da = float(gtex_da["FPR_at_0p5"])

tiles = [
    ("Toil accuracy @0.5", "0.515 -> 0.935", ADAPT, "deployed -> adapted"),
    ("Toil specificity", "0.03 -> 0.89", GOOD, "true-normal recovery"),
    ("Toil AUC", "0.992 -> 0.994", INK, "historical observation"),
    ("GTEx normals FPR", f"{gtex_fpr_dep:.3f} -> {gtex_fpr_da:.3f}", BAD, "all-normal: only partial"),
]

tile_html = "\n".join(
    f'<div class="tile"><div class="tl">{html.escape(t)}</div>'
    f'<div class="tv" style="color:{c}">{html.escape(v)}</div>'
    f'<div class="tc">{html.escape(cap)}</div></div>'
    for t, v, c, cap in tiles)

DOC = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Experimental cross-platform cohort standardization</title>
<style>
:root{{--ink:{INK};--mut:{MUT};--grid:{GRID};--bg:#ffffff;--card:#f7f9fc;--bd:#e3e8ef;}}
@media (prefers-color-scheme:dark){{:root{{--ink:#e7ecf3;--mut:#9aa5b1;--grid:#2a3340;--bg:#0f141b;--card:#161d27;--bd:#28313d;}}}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}}
.wrap{{max-width:940px;margin:0 auto;padding:34px 22px 60px}}
h1{{font-size:24px;margin:0 0 6px}} .lede{{color:var(--mut);margin:0 0 26px;font-size:15px}}
.tiles{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:0 0 30px}}
.tile{{background:var(--card);border:1px solid var(--bd);border-radius:12px;padding:14px 16px}}
.tl{{font-size:12px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}}
.tv{{font-size:22px;font-weight:700;margin:4px 0 2px;font-variant-numeric:tabular-nums}}
.tc{{font-size:12px;color:var(--mut)}}
.card{{background:var(--card);border:1px solid var(--bd);border-radius:14px;padding:16px 18px;margin:0 0 20px;overflow-x:auto}}
.ttl{{fill:var(--ink);font-size:14px;font-weight:700}} .sub{{fill:var(--mut);font-size:11px}}
.grid{{stroke:var(--grid);stroke-width:1}} .axl{{fill:var(--mut);font-size:10px}}
.cat{{fill:var(--ink);font-size:11px}} .val{{fill:var(--ink);font-size:10px;font-weight:600}}
.leg{{fill:var(--mut);font-size:11px}}
p{{margin:0 0 12px}} .note{{color:var(--mut);font-size:13px}}
code{{background:var(--card);border:1px solid var(--bd);border-radius:5px;padding:1px 5px;font-size:13px}}
.take{{border-left:3px solid {ADAPT};padding:6px 0 6px 14px;margin:16px 0}}
</style></head><body><div class="wrap">
<h1>Experimental cross-platform cohort standardization</h1>
<p class="lede">Historical benchmarks found severe threshold shift on foreign RNA-seq
pipelines. Standardizing on the input cohort is a transductive, composition-dependent
experiment &mdash; not probability calibration or a guarantee that threshold 0.5 transfers.</p>
<div class="tiles">{tile_html}</div>
<div class="card">{c1}</div>
<div class="take"><b>In the fixed historical Toil/RSEM cohort, accuracy changed from 0.515 to 0.935</b>
(specificity 0.03 &rarr; 0.89). This retrospective observation is not independent validation on a new batch.</div>
<div class="card">{c2}</div>
<p class="note">Cohort standardization recenters on the cohort mean, so it needs an internal tumor/normal
contrast. Results were stronger for mixed historical resamples and weak for
near-pure cohorts: on the all-normal GTEx panel it only moves FPR from {gtex_fpr_dep:.3f} to {gtex_fpr_da:.3f}
at 0.5. Do not adapt near-single-class cohorts. Labeled calibration is still apparent until independently confirmed.</p>
<p class="note">Reproduce: <code>python run_adaptation_benchmark.py</code> &nbsp;|&nbsp; Score:
<code>python cohort_adapt_score.py input.csv --adapt cohort_zscore</code> &nbsp;|&nbsp;
Frozen model reproduced to <code>max|&Delta;p| = 4.9e-7</code>. Generated 2026-07-05.</p>
</div></body></html>"""

out = os.path.join(HERE, "cross_platform_adaptation.html")
with open(
    out, "w", encoding="ascii", errors="xmlcharrefreplace", newline="\n"
) as f:
    f.write(DOC)
print("wrote", out, len(DOC), "bytes")
