# Bulk RNA-seq Analysis — NPM1 Primitive vs Committed (and CT45)

> **Scope:** This README documents the "BulkRNA" section of the Neoantigen_SingleCell project. It explains inputs, steps, outputs, and how to reproduce the results for bulk RNA-seq analyses (NPM1 primitive/committed scoring, CT45 associations, survival models).

---

## Table of contents

1. [Overview](#overview)
2. [Repository layout](#repository-layout)
3. [Data inputs](#data-inputs)
4. [Environment & dependencies](#environment--dependencies)
5. [Quick start](#quick-start)
6. [Methods](#methods)

   * [ID normalization (Ensembl → HGNC)](#id-normalization-ensembl--hgnc)
   * [Signature scoring (ssGSEA / weighted Z)](#signature-scoring-ssgsea--weighted-z)
   * [CT45 cluster association](#ct45-cluster-association)
   * [Survival analysis](#survival-analysis)
   * [Forest plot](#forest-plot)
   * [Multi‑panel figure assembly](#multi-panel-figure-assembly)
7. [Outputs](#outputs)
8. [Tips for reproducibility](#tips-for-reproducibility)
9. [Troubleshooting](#troubleshooting)
10. [Citations](#citations)

---

## Overview

This module annotates each bulk RNA‑seq sample as **NPM1‑primitive‑like** or **NPM1‑committed‑like**, relates those labels to **CT45 clusters**, and evaluates outcome differences with **Kaplan–Meier** and **Cox** models.

**Highlights**

* Robust signature scoring using **ssGSEA (GSVA)** with the current parameter‑object API.
* Safe handling of **gene identifiers** (Ensembl IDs → HGNC symbols; duplicate collapse).
* Clear association visuals (100% stacked bars) and **forest plots** for multivariable models.
* One‑click multi‑panel export via **patchwork**.

---

## Repository layout

A suggested structure (adapt names to your repo):

```
BulkRNA/
├─ README.md                 # this file
├─ BulkRNA.Rmd               # main narrative analysis (render to HTML/PDF)
├─ R/                        # helper functions used by the Rmd
│  ├─ ids.R                  # Ensembl→HGNC mapping & collapsing helpers
│  ├─ scoring.R              # ssGSEA / weighted Z scoring utilities
│  ├─ plots.R                # ggplot themes, forest/stacked-bar, KM wrappers
│  └─ utils.R                # small helpers (fmt_p, factor levels, etc.)
├─ data/
│  ├─ expression/            # input matrices (genes×samples; normalized)
│  ├─ meta/                  # clinical/metadata: age, sex, mutations, etc.
│  └─ signatures/            # meta-analysis table (GENENAME, estimate, fdr)
├─ results/
│  ├─ tables/
│  │  ├─ npm1_calls_*.csv
│  │  └─ cox_summary.csv
│  └─ figures/
│     ├─ assoc_ct45_vs_npm1.png
│     ├─ km_by_call.png
│     ├─ forest_cox.png
│     └─ all_panels.png
└─ renv/ (optional)          # if you use renv for reproducibility
```

> If you keep everything in a single `BulkRNA.Rmd`, that’s fine—consider moving long helper functions into `R/` to keep the Rmd narrative clean.

---

## Data inputs

* **Expression matrix**: normalized bulk RNA‑seq matrix (genes × samples). Gene rows may be **Ensembl** (with/without version) or **HGNC symbols**.
* **Meta table**: sample‑level annotations: `sample_id`, `NPM1`, `CT45 cluster` (e.g., `High_CT45`/`Low_CT45`), survival fields (`OS`, `Status`), covariates (`Age`, `Sex`, etc.), and driver mutations (e.g., `Mut_FLT3_ITD`).
* **Signature table**: meta‑analysis of NPM1 primitive vs committed with columns `GENENAME`, `estimate`, `fdr`.

> **Direction convention** used here: `estimate < 0` ⇒ **Primitive‑up**, `estimate > 0` ⇒ **Committed‑up**.

---

## Environment & dependencies

* **R ≥ 4.2** (tested with 4.3+)
* CRAN: `tidyverse`, `data.table`, `ggpubr`, `patchwork`, `survival`, `survminer`, `broom`
* Bioconductor: `GSVA (≥ 2.0 with new Param API)`, `org.Hs.eg.db`, `AnnotationDbi`

Optional: `renv` for project‑local libraries.

```r
# Core installs
install.packages(c("tidyverse","data.table","ggpubr","patchwork","survival","broom"))
BiocManager::install(c("GSVA","org.Hs.eg.db","AnnotationDbi","survminer"))
```

---

## Quick start

Render the analysis Rmd (adjust filename if needed):

```r
rmarkdown::render("BulkRNA/BulkRNA.Rmd")
```

Or source the key steps from the console:

```r
source("BulkRNA/R/ids.R")
source("BulkRNA/R/scoring.R")
source("BulkRNA/R/plots.R")
```

---

## Methods

### ID normalization (Ensembl → HGNC)

1. Strip Ensembl versions: `ENSG00000141510.18 → ENSG00000141510`.
2. Map to symbols with `org.Hs.eg.db`.
3. Collapse duplicate symbols by **mean** expression (bulk-friendly).

```r
ens <- sub("\\..*$","", rownames(expr))
map <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens), keytype="ENSEMBL", columns="SYMBOL")
map <- map[!is.na(map$SYMBOL), ]
expr_sym <- rowsum(expr[ens %in% map$ENSEMBL, , drop=FALSE], map$SYMBOL[match(ens, map$ENSEMBL)], reorder=FALSE)
expr_sym <- expr_sym / as.vector(table(map$SYMBOL))[rownames(expr_sym)]
rownames(expr_sym) <- toupper(rownames(expr_sym))
```

### Signature scoring (ssGSEA / weighted Z)

**Build sets from meta table** (flip signs to match convention):

```r
alpha <- 0.05
prim <- sig |> dplyr::filter(fdr <= alpha, estimate < 0) |> dplyr::pull(GENENAME) |> toupper()
comm <- sig |> dplyr::filter(fdr <= alpha, estimate > 0) |> dplyr::pull(GENENAME) |> toupper()
prim <- intersect(prim, rownames(expr_sym)); comm <- intersect(comm, rownames(expr_sym))
```

**ssGSEA (new GSVA API):**

```r
param <- GSVA::ssgseaParam(exprData = as.matrix(expr_sym),
                           geneSets = list(Primitive=prim, Committed=comm),
                           alpha = 0.25, normalize = TRUE)
ss <- GSVA::gsva(param)
score <- ss["Primitive", ] - ss["Committed", ]
call  <- ifelse(score > 0, "NPM1-primitive-like", "NPM1-committed-like")
```

**Weighted Z (uses effect sizes directly):**

```r
w <- setNames(sig$estimate, toupper(sig$GENENAME))
keep <- intersect(names(w), rownames(expr_sym))
Z <- t( scale(t(expr_sym[keep, , drop=FALSE])) )
score_z <- as.numeric(crossprod(matrix(w[keep], nrow=1), Z))
call_z  <- ifelse(score_z > 0, "NPM1-primitive-like", "NPM1-committed-like")
```

### CT45 cluster association

Compute 100% stacked bars and Fisher’s exact test:

```r
tab <- with(meta, table(cluster, call))
ft  <- fisher.test(tab)
assoc_df <- meta |>
  dplyr::filter(!is.na(call), !is.na(cluster)) |>
  dplyr::count(call, cluster) |>
  dplyr::group_by(call) |>
  dplyr::mutate(prop = n / sum(n), lab_y = cumsum(prop) - prop/2) |>
  dplyr::ungroup()

p_assoc <- ggplot(assoc_df, aes(call, prop, fill = cluster)) +
  geom_col(width = 0.75, color = "white") +
  geom_text(aes(y = lab_y, label = n), color = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(High_CT45 = "turquoise3", Low_CT45 = "indianred2")) +
  labs(y = "Proportion within NPM1 subgroup",
       subtitle = sprintf("Fisher’s p %s; OR = %.2f (95%% CI %.2f–%.2f)",
                          ifelse(ft$p.value < .001, "<0.001", sprintf("= %.3f", ft$p.value)),
                          unname(ft$estimate), ft$conf.int[1], ft$conf.int[2])) +
  theme_bw()
```

### Survival analysis

Kaplan–Meier by CT45 within each NPM1 call, plus Cox model.

```r
sf <- survfit(Surv(OS, Status) ~ cluster + call, data = meta)
km  <- survminer::ggsurvplot(sf, data = meta, facet.by = "call",
                             legend.title = "CT45 cluster",
                             palette = c("turquoise3","indianred2"),
                             risk.table = TRUE, risk.table.height = 0.18,
                             pval = TRUE)

fit <- coxph(Surv(OS, Status) ~ cluster + call + Age + Mut_FLT3_ITD + Mut_IDH2_p140 + Mut_IDH1 + Mut_TET2 + BlastsPB + LymphocytesPB + Sex,
             data = meta)
```

### Forest plot

```r
coefs <- summary(fit)$coefficients
ci    <- confint(fit)
forest_df <- data.frame(term = rownames(coefs),
                        HR = exp(coefs[,"coef"]),
                        lo = exp(ci[,1]), hi = exp(ci[,2]),
                        p  = coefs[,"Pr(>|z|)"])
forest_df$label <- sub("^cluster","CT45: ", forest_df$term)
forest_df$label[forest_df$term == "age"] <- "Age (per year)"
forest_df$label <- factor(forest_df$label, levels = rev(forest_df$label))

p_forest <- ggplot(forest_df, aes(HR, label)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_point() +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
  scale_x_log10(breaks = c(0.5,1,2,3,4)) +
  labs(x = "Hazard ratio (log scale)") + theme_bw()
```

### Multi‑panel figure assembly

Use **patchwork** and extract the survival curve + risk table.

```r
wrap_surv <- function(gsp, heights = c(3,1)) {
  if (!is.null(gsp$table)) patchwork::wrap_plots(list(gsp$plot, gsp$table), ncol = 1, heights = heights)
  else gsp$plot
}
sp_km <- wrap_surv(km)
final_fig <- ((p_assoc | p_forest) / sp_km) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")
```

---

## Outputs

* `results/tables/`

  * `npm1_calls_ssgsea.csv` or `npm1_calls_weightedZ.csv`
  * `cox_summary.csv`
* `results/figures/`

  * `assoc_ct45_vs_npm1.png`
  * `km_by_call.png` (with optional risk table)
  * `forest_cox.png`
  * `all_panels.png` (assembled)

---

## Tips for reproducibility

* Prefer **relative paths** (e.g., with `here::here()`), not absolute `C:/...`.
* Keep large raw data outside the repo; add sample/synthetic data for examples.
* Use `set.seed()` before resampling; capture `sessionInfo()` at the end of the Rmd.
* Consider `renv` for portable environments.

---

## Troubleshooting

* **GSVA error (defunct API):** Use `ssgseaParam(...)` then `gsva(param)`.
* **singscore error `argument is of length zero`:** One of the up/down sets is empty. Check overlaps with the matrix rownames.
* **Ensembl IDs not matching:** Strip versions and map to symbols before scoring.
* **`ggsurvplot` + patchwork errors:** Combine `$plot` and `$table` with `wrap_plots()` first.

---

## Citations

* Barbie DA et al., *Nature* 2009 — ssGSEA method.
* Hänzelmann S et al., *BMC Bioinformatics* 2013 — GSVA package.
* Therneau TM — `survival` and `coxph` documentation.
