#!/usr/bin/env python
"""Self-contained, JS-free HTML summary (inline SVG) for the cancer-type
classifier. ASCII-only so it renders identically from file://."""
import os
import csv
import json
import html

HERE = os.path.dirname(os.path.abspath(__file__))

def rows(name):
    with open(os.path.join(HERE, name), newline="") as f:
        return list(csv.reader(f))

cm_raw = rows("cancer_type_confusion_matrix.csv")
labels_file = cm_raw[0][1:]
cm = {r[0]: {labels_file[j]: int(v) for j, v in enumerate(r[1:])} for r in cm_raw[1:]}
per = {r["cancer_type"]: r for r in [dict(zip(rows("cancer_type_per_class_metrics.csv")[0], r))
       for r in rows("cancer_type_per_class_metrics.csv")[1:]]}
summ = json.load(open(os.path.join(HERE, "cancer_type_summary.json")))

# biologically-clustered order so confusions sit near the diagonal
ORDER = ["THCA","PRAD","BRCA","UCEC","BLCA","KIRC","KIRP","KICH","LIHC","CHOL",
         "LUAD","LUSC","HNSC","ESCA","STAD","COAD","READ"]

INK="#1f2933"; MUT="#8a94a6"; BD="#e3e8ef"
DIAG="#2f6df6"; ERR="#e8590c"

def lerp(c1, c2, t):
    return tuple(round(a + (b - a) * t) for a, b in zip(c1, c2))
def rgb(c): return f"rgb({c[0]},{c[1]},{c[2]})"

def heatmap():
    n = len(ORDER); cell = 27; padL = 52; padT = 60
    W = padL + n*cell + 14; H = padT + n*cell + 40
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="confusion matrix">']
    s.append(f'<text x="{padL}" y="18" class="ttl">Confusion matrix (patient-held-out, row-normalized)</text>')
    s.append(f'<text x="{padL}" y="34" class="sub">rows = true type, columns = predicted; number = sample count. Types ordered to group related tissues.</text>')
    for i, tr in enumerate(ORDER):
        rowtot = sum(cm[tr].values()) or 1
        y = padT + i*cell
        s.append(f'<text x="{padL-6}" y="{y+cell/2+3:.0f}" class="hlab" text-anchor="end">{tr}</text>')
        s.append(f'<text x="{padL+i*cell+cell/2:.0f}" y="{padT-6}" class="hlab" text-anchor="middle" transform="rotate(-90 {padL+i*cell+cell/2:.0f} {padT-8})">{tr}</text>')
        for j, tc in enumerate(ORDER):
            v = cm[tr][tc]; frac = v / rowtot; x = padL + j*cell
            if v == 0:
                fill = "#f4f7fb"
            elif i == j:
                fill = rgb(lerp((225,235,252), (47,109,246), min(1, frac)))
            else:
                fill = rgb(lerp((253,240,230), (232,89,12), min(1, frac*1.6)))
            s.append(f'<rect x="{x}" y="{y}" width="{cell-1}" height="{cell-1}" fill="{fill}"/>')
            if v:
                tcol = "#fff" if (i==j and frac>0.55) or (i!=j and frac*1.6>0.55) else INK
                s.append(f'<text x="{x+cell/2-0.5:.0f}" y="{y+cell/2+3:.0f}" class="cnum" fill="{tcol}" text-anchor="middle">{v}</text>')
    s.append('</svg>')
    return "".join(s)

def f1bars():
    items = sorted(ORDER, key=lambda t: float(per[t]["f1"]), reverse=True)
    W=560; rowh=20; padL=54; padT=44; H=padT+len(items)*rowh+10
    barW=W-padL-70
    s=[f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="per-type F1">']
    s.append(f'<text x="{padL}" y="18" class="ttl">Per-type F1 (patient-held-out)</text>')
    s.append(f'<text x="{padL}" y="33" class="sub">READ collapses into COAD (colorectal = one disease); unique-marker tissues are perfect.</text>')
    for i,t in enumerate(items):
        f1=float(per[t]["f1"]); y=padT+i*rowh
        col=rgb(lerp((232,89,12),(47,109,246), f1))  # low=orange, high=blue
        s.append(f'<text x="{padL-6}" y="{y+rowh/2+3:.0f}" class="hlab" text-anchor="end">{t}</text>')
        s.append(f'<rect x="{padL}" y="{y+2:.0f}" width="{barW}" height="{rowh-6}" fill="#eef2f7"/>')
        s.append(f'<rect x="{padL}" y="{y+2:.0f}" width="{barW*f1:.0f}" height="{rowh-6}" rx="2" fill="{col}"/>')
        s.append(f'<text x="{padL+barW+6}" y="{y+rowh/2+3:.0f}" class="cnum" fill="{INK}" text-anchor="start">{f1:.2f} (n={per[t]["support"]})</text>')
    s.append('</svg>')
    return "".join(s)

tiles=[("Accuracy", f'{summ["accuracy"]:.3f}', DIAG, "patient-held-out"),
       ("Balanced accuracy", f'{summ["balanced_accuracy"]:.3f}', INK, "mean over 17 types"),
       ("Macro F1", f'{summ["macro_f1"]:.3f}', INK, "unweighted"),
       ("Perfect types", "THCA, PRAD", "#12b886", "F1 = 1.00")]
tile_html="".join(f'<div class="tile"><div class="tl">{html.escape(t)}</div>'
    f'<div class="tv" style="color:{c}">{html.escape(v)}</div><div class="tc">{html.escape(cap)}</div></div>'
    for t,v,c,cap in tiles)

DOC=f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TCGA cancer-type classifier</title><style>
:root{{--ink:{INK};--mut:{MUT};--bd:{BD};--bg:#fff;--card:#f7f9fc;}}
@media (prefers-color-scheme:dark){{:root{{--ink:#e7ecf3;--mut:#9aa5b1;--bd:#28313d;--bg:#0f141b;--card:#161d27;}}}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}}
.wrap{{max-width:940px;margin:0 auto;padding:34px 22px 60px}}
h1{{font-size:23px;margin:0 0 6px}} .lede{{color:var(--mut);margin:0 0 24px}}
.tiles{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:0 0 26px}}
.tile{{background:var(--card);border:1px solid var(--bd);border-radius:12px;padding:13px 15px}}
.tl{{font-size:12px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}}
.tv{{font-size:21px;font-weight:700;margin:3px 0 2px;font-variant-numeric:tabular-nums}}
.tc{{font-size:12px;color:var(--mut)}}
.card{{background:var(--card);border:1px solid var(--bd);border-radius:14px;padding:16px 18px;margin:0 0 20px;overflow-x:auto}}
.ttl{{fill:var(--ink);font-size:13px;font-weight:700}} .sub{{fill:var(--mut);font-size:10.5px}}
.hlab{{fill:var(--ink);font-size:10px}} .cnum{{font-size:10px;font-weight:600}}
.take{{border-left:3px solid {ERR};padding:6px 0 6px 14px;margin:14px 0;color:var(--ink)}}
.note{{color:var(--mut);font-size:13px}} code{{background:var(--card);border:1px solid var(--bd);border-radius:5px;padding:1px 5px;font-size:13px}}
</style></head><body><div class="wrap">
<h1>TCGA cancer-type (tissue-of-origin) classifier</h1>
<p class="lede">Multinomial logistic regression over 1,000 genes classifies a tumor's bulk RNA-seq
profile into 17 TCGA cancer types. 1,440 tumors, 1,438 patients; evaluated patient-held-out.</p>
<div class="tiles">{tile_html}</div>
<div class="card">{heatmap()}</div>
<div class="take"><b>The errors are biologically adjacent tissues, not noise.</b> The bright off-diagonal
block is READ&rarr;COAD (15 of 20): rectal and colon cancer are one disease. Other confusions cluster
within kidney (KIRC/KIRP/KICH), lung/squamous (LUAD/LUSC/HNSC), upper-GI (ESCA/STAD) and hepatobiliary (LIHC/CHOL).</div>
<div class="card">{f1bars()}</div>
<p class="note">Marker genes recover canonical tissue markers (THCA: TG/TPO/FOXE1; PRAD: ACP3/NKX3-1/KLK4;
LIHC: SHBG/AMBP/GC; LUSC: SFTPB/SFTPA1), confirming genuine tissue-of-origin signal.
Score new samples: <code>python predict_cancer_type.py input.csv</code>. Pure-numpy model reproduces
scikit-learn to max|&Delta;p|=1.5e-8. Generated 2026-07-06.</p>
</div></body></html>"""

out=os.path.join(HERE,"cancer_type_classifier.html")
open(out,"w",encoding="ascii",errors="xmlcharrefreplace").write(DOC)
print("wrote", out, len(DOC), "bytes")
