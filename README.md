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

1. **QC**: Quality check the reads
2. **Merge**: Merge them as single end reads
3. **Right‑size select**: Additional QC to check we have right size amplicons
4. **Importing the reads** :Create a tab delimited file which has filepath for reads
5. **dereplicate**: and do **de novo 99% clustering**.
6. **References**: download Norovirus to create database, re‑header (add G/P labels).
7. **Train classifier** and **test classifier** on sample reps.
8. **Assign taxonomy** (G‑type, P‑type).
9. **Relative abundance**: executed notebook summarizes genotype/P‑type composition across samples.

---
---

## 📈 Pipeline at‑a‑glance (arrows)

```
Sample Screening and Library Preparation
Wastewater → RT‑dPCR screening → 2‑step semi‑nested PCR (BC junction) → Illumina 2×300 → FASTQs

Downstream Analysis
R1/R2 FASTQs
  ↓
Quality Control
  ↓
Merge  
  ↓
import reads 
  ↓
VSEARCH dereplicate → VSEARCH de‑novo cluster @99%
  ↓
Visualize (table , representative seqs for each cluster )
  ↓
refs: download → re‑header (add G/P labels) → extract BC (primer sets) → merge refs
  ↓
Train classifier (Naive Bayes) → classify‑sklearn
  ↓
Taxonomy assigned (G‑type, P‑type)
  ↓
Export: feature‑table → feature‑table‑frequency  +  taxonomy 
  ↓
Abundance calculation using python → (GI/GII; G/P types)




```

---


## What we observed 

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
