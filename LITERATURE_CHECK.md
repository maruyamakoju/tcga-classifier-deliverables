# Literature check — do the top genes match published pan-cancer signatures?

**Verdict: yes, strongly — and the picture is richer than "stromal/ECM remodeling."**
The classifier separates tumor from normal using a *bidirectional* signature: it keys on
genes that are **gained in tumors** (ECM remodeling / cancer-associated-fibroblast /
desmoplasia programs) **and** on normal-tissue architecture / microvascular / tumor-
suppressor genes that are **lost in tumors**. Both arms are well documented pan-cancer,
which is exactly why the model generalizes to cancer types it never trained on (this
gain/loss pattern is shared across epithelial cancers).

## Arm 1 — gained in tumors (top drivers of the "tumor" call)

- **COL10A1** (among the top logistic-regression coefficients, ranked #3 after
  BMPER and GABRD): a pan-cancer hallmark. Reported
  significantly overexpressed in **13 cancer types** — bladder, breast, cholangio, colon,
  esophageal, head & neck, LUAD, LUSC, prostate, rectal, stomach, thyroid, endometrial —
  essentially the same set used here. ECM structural collagen that promotes invasion/metastasis.
- **CTHRC1**: identified as a **major pan-cancer ECM (matrisome) regulator** in a dedicated
  pan-cancer analysis across TCGA.
- **COMP** and **SPP1** (osteopontin): both members of a validated **9-gene pan-cancer ECM
  signature** associated with poor survival across lung, gastric, ovarian and colon cancers.
- **MMP11** (top XGBoost gene) and **MMP9**: matrix metalloproteinases; MMP11 is
  cancer-associated-fibroblast–derived and upregulated pan-cancer, driving ECM breakdown,
  invasion and metastasis.
- **PDGFRA** (top XGBoost gene): a classic **cancer-associated-fibroblast** marker,
  upregulated in cancer stroma. Caveat: at single-cell resolution PDGFRA marks only a
  *subset* of CAFs and is also expressed by normal fibroblasts/pericytes, so at bulk level
  its tumor-vs-normal direction is tissue-dependent — treat it as a stromal-activation
  proxy, not a clean tumor marker.

## Arm 2 — lost in tumors (top drivers of the "normal" call)

- **DCN** (decorin): a well-established **tumor suppressor**; TCGA shows significant
  reduction across **18 solid tumors**. Restrains EGFR/Met/IGF-1R/VEGFR2 signaling and
  lymphangiogenesis; its loss activates pro-tumor stroma.
- **LYVE1**: lymphatic-endothelial marker, significantly **lower in tumor** than
  para-cancerous tissue — loss of normal lymphatic architecture.
- **TGFBR3** (TGF-β type III receptor): reduced in breast, kidney, lung, ovary, prostate
  and liver cancers; a suppressor of invasion, proliferation and angiogenesis.
- **FBLN5** (fibulin-5): ~5-fold downregulated in ovarian cancer; suppresses MMP9,
  angiogenesis and cell motility.
- **MMRN1**, **ABCA8**, **DPT** (dermatopontin), **TNXB**, **MFAP4**, **ANGPTL1**: normal
  stromal-matrix and microvascular markers, consistent with the same "loss of normal
  tissue architecture" theme (MMRN family are endothelial/vascular markers).

## Interpretation

The original report's claim — "top genes dominated by stromal/ECM remodeling markers,
consistent with a known pan-cancer signature" — **holds up**, and the more precise reading
is: **tumor = (desmoplastic ECM / CAF program UP) + (normal tissue, lymphatic/vascular and
tumor-suppressor markers DOWN).** This bidirectional, tissue-agnostic contrast is the
mechanistic explanation for the leave-one-cancer-type-out result (AUC 0.994): the model
isn't memorizing per-cancer patterns, it's reading a shared malignant-vs-normal tissue
state.

## Caveat carried over

TCGA "normal" is tumor-adjacent resected tissue, so part of Arm 2 reflects **field/adjacent
normal** rather than healthy-population tissue. The signature is therefore validated as
tumor-vs-adjacent-normal; generalization to healthy-donor tissue is untested.

## Sources
- COL10A1 pan-cancer review: https://pmc.ncbi.nlm.nih.gov/articles/PMC11487528/ ; LUAD ECM remodeling: https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2020.573534/full
- CTHRC1 pan-cancer matrisome: https://pmc.ncbi.nlm.nih.gov/articles/PMC9529084/ (PubMed https://pubmed.ncbi.nlm.nih.gov/36190948/)
- Matrisome / tumor-microenvironment review (SPP1, COMP signature context): https://www.sciencedirect.com/science/article/pii/S0304419X24001094
- MMP11 CAF-derived, pan-cancer: https://pmc.ncbi.nlm.nih.gov/articles/PMC11875780/ ; MMP11/SPP1 biomarkers: https://pmc.ncbi.nlm.nih.gov/articles/PMC12082244/
- PDGFRA as CAF marker (pan-cancer): https://pmc.ncbi.nlm.nih.gov/articles/PMC9402225/ ; CAF heterogeneity: https://pmc.ncbi.nlm.nih.gov/articles/PMC9327514/
- Decorin (DCN) tumor suppressor, reduced in 18 solid tumors: https://www.nature.com/articles/s42003-020-01590-0 ; lymphangiogenesis: https://www.pnas.org/doi/10.1073/pnas.2317760121
- LYVE1 lower in tumor: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5698955/
- FBLN5 downregulated: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5865667/
- TGFBR3 tumor suppressor: https://www.mdpi.com/2072-6694/12/6/1375 ; HCC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8037431/
