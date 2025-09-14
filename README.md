# noro\_ww\_pipeline — Wastewater Surveillance of Norovirus

**One‑liner:** A short‑read amplicon pipeline to type **Norovirus** in **wastewater** by sequencing the **RdRp→VP1 (BC) junction**, enabling simultaneous **G‑type (VP1)** and **P‑type (RdRp)** calls.

---

## Why Norovirus & why wastewater?

* **Norovirus matters.** Norovirus is the major cause of gastroenteritis cases worldwide and in the USA. Globally, Norovirus is responsible for approximately 18% of all clinical cases of gastroenteritis and an estimated 207,000 deaths each year. 
* **Wastewater works.** Wastewater-based epidemiology combined with genomic sequencing, allows the tracking viral strains such SARS-CoV-2, Polio, RSV etc. in real time, allowing for early detection and better public health response. However, there remains a significant gap in applying wastewater surveillance methods for Norovirus. 

---

## Why the **BC (RdRp→VP1) junction**?

Norovirus genomes encode non‑structural proteins in **ORF1** (including **RdRp**) and capsid proteins in **ORF2** (**VP1**). The **BC junction** lies at the interface of RdRp and VP1 and is a hotspot for recombination. By amplifying this junction (\~**579 nt** in **GI**; \~**570 nt** in **GII**), a **single merged read** carries both:

* **P‑type** information (from **RdRp**, ORF1) and
* **G‑type** information (from **VP1**, ORF2).

This pairing enables dual typing of Norovirus.




## Our innovation

**Wet‑lab input to this pipeline**

* **dPCR screening** of wastewater extracts.
* **PCR** targeting the BC region (GI/GII assays).
* **Illumina 2×300 bp** sequencing to allow full amplicon recovery.

**Dry‑lab pipeline (this repo)**

1. **QC**: fastp → Trim Galore (Q20, error=0.2; Nextera/Illumina adapters).
2. **Merge**: BBMerge; record insert‑size histogram.
3. **Right‑size select**: keep **570–580 nt** merged reads (junction amplicon).
4. **QIIME 2 / VSEARCH path**: import (single‑end), **dereplicate**, **de novo 99% clustering**.
5. **References**: download HUCAT, re‑header (add G/P labels), **extract BC region** by multiple primer sets, merge refs.
6. **Train classifier** (Naive Bayes) on BC‑region refs + taxonomy; **test classifier** on sample reps.
7. **Assign taxonomy** (G‑type, P‑type).
8. **Relative abundance**: executed notebook summarizes genotype/P‑type composition across samples.

---
---

## 📈 Pipeline at‑a‑glance (arrows)

```
WET LAB
Wastewater → RT‑dPCR screening → 2‑step semi‑nested PCR (BC junction) → Illumina 2×300 → FASTQs

DRY LAB 
R1/R2 FASTQs
  ↓
fastp 
  ↓
Trim Galore 
  ↓
BBMerge 
  ↓
reformat.sh → keep 570–580 nt  
  ↓
QIIME 2 import (single‑end manifest)
  ↓
VSEARCH dereplicate → VSEARCH de‑novo cluster @99%
  ↓
Visualize (table.qzv, rep‑seqs.qzv)
  ↓
refs: download → re‑header (add G/P labels) → extract BC (primer sets) → merge refs
  ↓
Train classifier (Naive Bayes) → classify‑sklearn
  ↓
Taxonomy assigned (G‑type, P‑type)
  ↓
Export: feature‑table.biom → feature‑table‑frequency.tsv  +  taxonomy.tsv
  ↓
Notebook execution → (GI/GII; G/P types)




```

---


## What we observed (pilot summary)

* In **Genogroup I (GI)** samples, the most abundant variant was **GI.1 P\[1]**.
* In **Genogroup II (GII)** samples, the most abundant variant was **GII.4 P\[39]**.

---

## Long‑term vision

* **Routine genotype‑level WBE:** near‑real‑time dashboards of G/P types across catchments.
* **Benchmarking:** side‑by‑side ASV path for higher resolution and validation.

---

## What’s in this repo (minimal)

* **Snakemake workflow** from FASTQs → merged junction reads → QIIME 2/VSEARCH typing.
* **Env files** for reproducible installs; **config** with auto‑discovery of multiple samples.
* **Ref utilities** to prepare HUCAT‑based BC‑region training data.
* **Notebook** to compute **relative abundance** of G/P types.




**Inputs:** paired‑end FASTQs in `data/` (auto‑discovered or listed in `config.yaml`).
**Outputs:** merged 570–580 nt reads, QIIME artifacts/visualizations, feature table + taxonomy TSVs, executed notebook.



**Contact:** Anirudh Bhatia · [bhatiaan@oregonstate.edu](mailto:bhatiaan@oregonstate.edu)
**License:** MIT (code). Reference data may carry their own terms.
