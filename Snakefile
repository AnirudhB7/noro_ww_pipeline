import os
from pathlib import Path
from snakemake.shell import shell
shell.executable("bash")
shell.prefix("set -euo pipefail; ")

configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())
THREADS = config.get("threads", 8)
ADAPTERS = config.get("adapters", "illumina")  # or "nextera"
TG_E = config.get("trimgalore_error", 0.2)
TG_Q = config.get("trimgalore_quality", 20)
MINLEN = config.get("min_length", 50)
MINOV = config.get("bbmerge_minoverlap", 12)
LEN_MIN = config.get("len_select_min", 570)
LEN_MAX = config.get("len_select_max", 580)

# Reference + classifier config
HUCAT_FASTA = config["refs"]["hucat_fasta"]
NEW_HEADERS = config["refs"]["new_headers_mapping"]
UPDATED_FASTA = "refs/Updated_HUCAT_TAX_ACC_CHANGED.fasta"
TAXON_TSV = config["refs"]["taxonomy_tsv"]
PRIMER_SETS = config["refs"]["primer_sets_yaml"] # YAML file listing primer pairs

# Helper for Trim Galore adapter flag
TG_FLAG = "--nextera" if ADAPTERS.lower() == "nextera" else "--illumina"

rule all:
    input:
        # Core per-sample outputs through clustering
        expand("results/bbmerge/{s}/merged_{s}.fastq.gz", s=SAMPLES),
        expand("results/bbmerge/{s}/merge_stats_{s}.txt", s=SAMPLES),
        expand("results/merged_reads/{s}.fastq.gz", s=SAMPLES),
        "results/qiime/import/single-end-demux.qza",
        "results/qiime/derep_cluster/table-dn-99.qza",
        "results/qiime/derep_cluster/rep-seqs-dn-99.qza",
        "results/qiime/viz/table-dn-99.qzv",
        "results/qiime/viz/rep-seqs-dn-99.qzv",
        # Reference artifacts + classifier + taxonomy
        "results/qiime/refs/custom_capsid_db.qza",
        "results/qiime/refs/ref-seqs-merged.qza",
        "results/qiime/refs/custom_ref-taxonomy.qza",
        "results/qiime/refs/classifier.qza",
        "results/qiime/derep_cluster/taxonomy.qza",
        "results/qiime/viz/taxonomy.qzv",
        # Data Exports
        "results/qiime/export/exported_table/feature-table.biom",
        "results/qiime/export/taxonomy/taxonomy.tsv",
        "results/qiime/export/feature-table-frequency.tsv",
        # QC summary
        "results/qc/QC_SUMMARY.tsv"
        #Analysis
        "results/analysis/relative_abundance.executed.ipynb"
# ----------
# fastp
# ----------
rule fastp:
    conda: "envs/fastp.yaml"
    threads: THREADS
    input:
        R1=lambda w: config["samples"][w.s]["R1"],
        R2=lambda w: config["samples"][w.s]["R2"]
    output:
        R1="results/fastp/{s}/R1.fastq.gz",
        R2="results/fastp/{s}/R2.fastq.gz",
        html="results/fastp/{s}/fastp_report.html",
        json="results/fastp/{s}/fastp_report.json"
    params:
        q=TG_Q,
        minlen=MINLEN
    shell:
        r"""
        mkdir -p results/fastp/{wildcards.s}
        fastp \
          -i {input.R1} -I {input.R2} \
          -o {output.R1} -O {output.R2} \
          --detect_adapter_for_pe \
          --cut_front --cut_tail \
          --qualified_quality_phred {params.q} \
          --length_required {params.minlen} \
          -g \
          --html {output.html} \
          --json {output.json} \
          --thread {threads}
        """

# -------------------
# Trim Galore
# -------------------
rule trimgalore:
    conda: "envs/trimgalore.yaml"
    threads: THREADS
    input:
        R1="results/fastp/{s}/R1.fastq.gz",
        R2="results/fastp/{s}/R2.fastq.gz"
    output:
        R1="results/trimgalore/{s}/R1_val_1.fq.gz",
        R2="results/trimgalore/{s}/R2_val_2.fq.gz",
        r1_html="results/trimgalore/{s}/R1_val_1_fastqc.html",
        r2_html="results/trimgalore/{s}/R2_val_2_fastqc.html",
        r1_rep="results/trimgalore/{s}/R1_trimming_report.txt",
        r2_rep="results/trimgalore/{s}/R2_trimming_report.txt"
    params:
        flag=TG_FLAG,
        e=TG_E,
        q=TG_Q
    shell:
        r"""
        mkdir -p results/trimgalore/{wildcards.s}
        trim_galore {params.flag} -e {params.e} --quality {params.q} \
          --paired --fastqc --output_dir results/trimgalore/{wildcards.s} \
          {input.R1} {input.R2}

        # Normalize names (symlink to canonical outputs)
        R1_OUT=$(ls results/trimgalore/{wildcards.s}/*_val_1.fq.gz | head -n1)
        R2_OUT=$(ls results/trimgalore/{wildcards.s}/*_val_2.fq.gz | head -n1)
        R1_HTML=$(ls results/trimgalore/{wildcards.s}/*_val_1_fastqc.html | head -n1)
        R2_HTML=$(ls results/trimgalore/{wildcards.s}/*_val_2_fastqc.html | head -n1)
        R1_REP=$(ls results/trimgalore/{wildcards.s}/*R1*report.txt 2>/dev/null || true)
        R2_REP=$(ls results/trimgalore/{wildcards.s}/*R2*report.txt 2>/dev/null || true)
        ln -sf "$R1_OUT" {output.R1}
        ln -sf "$R2_OUT" {output.R2}
        ln -sf "$R1_HTML" {output.r1_html}
        ln -sf "$R2_HTML" {output.r2_html}
        
        [ -n "$R1_REP" ] && ln -sf "$R1_REP" {output.r1_rep} || echo "NA" > {output.r1_rep}
        [ -n "$R2_REP" ] && ln -sf "$R2_REP" {output.r2_rep} || echo "NA" > {output.r2_rep}
        """

# --------------
# BBMerge
# --------------
rule bbmerge:
    conda: "envs/bbtools.yaml"
    threads: THREADS
    input:
        R1="results/trimgalore/{s}/R1_val_1.fq.gz",
        R2="results/trimgalore/{s}/R2_val_2.fq.gz"
    output:
        merged="results/bbmerge/{s}/merged_{s}.fastq.gz",
        u1="results/bbmerge/{s}/unmerged_R1_{s}.fastq.gz",
        u2="results/bbmerge/{s}/unmerged_R2_{s}.fastq.gz",
        ihist="results/bbmerge/{s}/merge_stats_{s}.txt"
    params:
        minov=MINOV
    shell:
        r"""
        mkdir -p results/bbmerge/{wildcards.s}
        bbmerge.sh in1={input.R1} in2={input.R2} \
          out={output.merged} outu1={output.u1} outu2={output.u2} \
          ihist={output.ihist} minoverlap={params.minov} threads={threads}
        """

# -----------------------------
# Length-select merged
# -----------------------------
rule length_select_merged:
    conda: "envs/bbtools.yaml"
    input:
        merged="results/bbmerge/{s}/merged_{s}.fastq.gz"
    output:
        out="results/merged_reads/{s}.fastq.gz"
    params:
        lmin=LEN_MIN,
        lmax=LEN_MAX
    shell:
        r"""
        mkdir -p results/merged_reads
        reformat.sh in={input.merged} out={output.out} minlength={params.lmin} maxlength={params.lmax} ow=t
        """

# ---------------------------------
# QIIME2 import (single-end)
# ---------------------------------
rule make_manifest:
    output:
        manifest="results/qiime/import/single_end_phred33_manifest.tsv"
    input:
        expand("results/merged_reads/{s}.fastq.gz", s=SAMPLES)
    run:
        import os
        os.makedirs("results/qiime/import", exist_ok=True)
        with open(output.manifest, "w") as f:
            f.write("sample-id\tabsolute-filepath\n")
            for s in SAMPLES:
                abspath = os.path.abspath(f"results/merged_reads/{s}.fastq.gz")
                f.write(f"{s}\t{abspath}\n")

rule qiime_import_se:
    conda: ""envs/qiime2-amplicon-2024.10.yml""
    input:
        manifest="results/qiime/import/single_end_phred33_manifest.tsv"
    output:
        qza="results/qiime/import/single-end-demux.qza"
    shell:
        r"""
        qiime tools import \
          --type 'SampleData[SequencesWithQuality]' \
          --input-path {input.manifest} \
          --output-path {output.qza} \
          --input-format SingleEndFastqManifestPhred33V2
        """

# ------------------------------
# Dereplicate + cluster
# ------------------------------
rule dereplicate:
    conda: ""envs/qiime2-amplicon-2024.10.yml""
    input:
        qza="results/qiime/import/single-end-demux.qza"
    output:
        table="results/qiime/derep_cluster/table.qza",
        reps="results/qiime/derep_cluster/rep-seqs.qza"
    shell:
        r"""
        mkdir -p results/qiime/derep_cluster
        qiime vsearch dereplicate-sequences \
          --i-sequences {input.qza} \
          --o-dereplicated-table {output.table} \
          --o-dereplicated-sequences {output.reps}
        """

rule cluster_99:
    conda: ""envs/qiime2-amplicon-2024.10.yml""
    input:
        table="results/qiime/derep_cluster/table.qza",
        reps="results/qiime/derep_cluster/rep-seqs.qza"
    output:
        table99="results/qiime/derep_cluster/table-dn-99.qza",
        reps99="results/qiime/derep_cluster/rep-seqs-dn-99.qza"
    shell:
        r"""
        qiime vsearch cluster-features-de-novo \
          --i-table {input.table} \
          --i-sequences {input.reps} \
          --p-perc-identity 0.99 \
          --o-clustered-table {output.table99} \
          --o-clustered-sequences {output.reps99}
        """

# ------------------------
# Visualizations
# ------------------------
rule viz_feature_table:
    conda: ""envs/qiime2-amplicon-2024.10.yml""
    input:
        table="results/qiime/derep_cluster/table-dn-99.qza"
    output:
        qzv="results/qiime/viz/table-dn-99.qzv"
    shell:
        r"""
        mkdir -p results/qiime/viz
        qiime feature-table summarize \
          --i-table {input.table} \
          --o-visualization {output.qzv}
        """

rule viz_rep_seqs:
    conda: ""envs/qiime2-amplicon-2024.10.yml""
    input:
        reps="results/qiime/derep_cluster/rep-seqs-dn-99.qza"
    output:
        qzv="results/qiime/viz/rep-seqs-dn-99.qzv"
    shell:
        r"""
        qiime feature-table tabulate-seqs \
          --i-data {input.reps} \
          --o-visualization {output.qzv}
        """

# --------------------------------------
# HUCAT refs to reheader to import to extract BC
# --------------------------------------
rule reheader_fasta:
    input:
        fasta=HUCAT_FASTA,
        mapping=NEW_HEADERS
    output:
        fasta=UPDATED_FASTA
    run:
        # mapping file must have one header per line, including leading '>'
        mapping = {}
        with open(input.mapping) as f:
            for line in f:
                line=line.strip()
                if not line: continue
                hdr = line[1:] if line.startswith('>') else line
                key = hdr.split('_')[0]
                mapping[key]=hdr
        with open(input.fasta) as fin, open(output.fasta,'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    acc = line[1:].split()[0]
                    if acc in mapping:
                        fout.write('>' + mapping[acc] + '\n')
                    else:
                        fout.write(line)
                else:
                    fout.write(line)

rule import_ref_seqs:
    conda: ""envs/qiime2-amplicon-2024.10.yml""
    input:
        fasta=UPDATED_FASTA
    output:
        qza="results/qiime/refs/custom_capsid_db.qza"
    shell:
        r"""
        mkdir -p results/qiime/refs
        qiime tools import \
          --type 'FeatureData[Sequence]' \
          --input-path {input.fasta} \
          --output-path {output.qza}
        """

# Primer sets for extract-reads are defined in refs/primer_sets.yaml
#Norovirus round 2 primers for extracting BC region of Norovirus



import yaml
with open(PRIMER_SETS) as fh:
    PRIMERS = yaml.safe_load(fh)["sets"]

rule extract_reads:
    conda: envs/qiime2-amplicon-2024.10.yml"
    input:
        db="results/qiime/refs/custom_capsid_db.qza"
    output:
        qza="results/qiime/refs/{name}.qza"
    params:
        fwd=lambda w: next(p['fwd'] for p in PRIMERS if p['name']==w.name),
        rev=lambda w: next(p['rev'] for p in PRIMERS if p['name']==w.name)
    shell:
        r"""
        qiime feature-classifier extract-reads \
          --i-sequences {input.db} \
          --p-f-primer {params.fwd} \
          --p-r-primer {params.rev} \
          --o-reads {output.qza}
        """
rule merge_ref_seqs:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        expand("results/qiime/refs/{name}.qza", name=[p['name'] for p in PRIMERS])
    output:
        merged="results/qiime/refs/ref-seqs-merged.qza"
    run:
        args = " ".join(f"--i-data {p}" for p in input)
        shell(f"""
            qiime feature-table merge-seqs \
              {args} \
              --o-merged-data {output.merged}
        """)


# -------------------------
# Train classifier
# -------------------------
rule import_taxonomy:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        tsv=TAXON_TSV
    output:
        tax="results/qiime/refs/custom_ref-taxonomy.qza"
    shell:
        r"""
        qiime tools import \
          --type 'FeatureData[Taxonomy]' \
          --input-format HeaderlessTSVTaxonomyFormat \
          --input-path {input.tsv} \
          --output-path {output.tax}
        """

rule train_nb_classifier:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        seqs="results/qiime/refs/ref-seqs-merged.qza",
        tax="results/qiime/refs/custom_ref-taxonomy.qza"
    output:
        clf="results/qiime/refs/classifier.qza"
    shell:
        r"""
        qiime feature-classifier fit-classifier-naive-bayes \
          --i-reference-reads {input.seqs} \
          --i-reference-taxonomy {input.tax} \
          --o-classifier {output.clf}
        """

# ---------------------
# Classify OTUs
# ---------------------
rule classify_reads:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        clf="results/qiime/refs/classifier.qza",
        reps="results/qiime/derep_cluster/rep-seqs-dn-99.qza"
    output:
        tax="results/qiime/derep_cluster/taxonomy.qza"
    shell:
        r"""
        qiime feature-classifier classify-sklearn \
          --i-classifier {input.clf} \
          --i-reads {input.reps} \
          --o-classification {output.tax}
        """

rule viz_taxonomy:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        tax="results/qiime/derep_cluster/taxonomy.qza"
    output:
        qzv="results/qiime/viz/taxonomy.qzv"
    shell:
        r"""
        qiime metadata tabulate \
          --m-input-file {input.tax} \
          --o-visualization {output.qzv}
        """

# -------------------
# Data Exports
# -------------------
rule export_table:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        table="results/qiime/derep_cluster/table-dn-99.qza",
        tax="results/qiime/derep_cluster/taxonomy.qza"
    output:
        biom="results/qiime/export/exported_table/feature-table.biom",
        tax_tsv="results/qiime/export/taxonomy/taxonomy.tsv"
    shell:
        r"""
        mkdir -p results/qiime/export/exported_table results/qiime/export/taxonomy
        qiime tools export --input-path {input.table} --output-path results/qiime/export/exported_table
        qiime tools export --input-path {input.tax}   --output-path results/qiime/export/taxonomy
        """

rule biom_to_tsv:
    conda: "envs/qiime2-amplicon-2024.10.yml"
    input:
        biom="results/qiime/export/exported_table/feature-table.biom"
    output:
        tsv="results/qiime/export/feature-table-frequency.tsv"
    shell:
        r"""
        biom convert -i {input.biom} -o {output.tsv} --to-tsv
        """

# -------------------
# QC summary
# -------------------
rule qc_summary:
    input:
        expand("results/fastp/{s}/fastp_report.json", s=SAMPLES),
        expand("results/trimgalore/{s}/R1_trimming_report.txt", s=SAMPLES),
        expand("results/trimgalore/{s}/R2_trimming_report.txt", s=SAMPLES),
        expand("results/bbmerge/{s}/merge_stats_{s}.txt", s=SAMPLES)
    output:
        "results/qc/QC_SUMMARY.tsv"
    run:
        import json, re, csv, os
        os.makedirs("results/qc", exist_ok=True)
        def parse_tg(path):
            tot=kept=None
            if not os.path.exists(path):
                return (None,None)
            with open(path) as fh:
                for ln in fh:
                    if "Total reads processed:" in ln:
                        tot=int(ln.strip().split()[-1].replace(',',''))
                    if "Reads written (passing filters):" in ln:
                        kept=int(ln.strip().split()[-1].replace(',',''))
            return (tot,kept)
        rows=[["sample","fastp_in","fastp_out","tg_r1_in","tg_r1_kept","tg_r2_in","tg_r2_kept","bbmerge_pct_merged","insert_mean","insert_median","insert_mode"]]
        for s in SAMPLES:
            with open(f"results/fastp/{s}/fastp_report.json") as f:
                js=json.load(f)
            fp_in=js["summary"]["before_filtering"]["total_reads"]
            fp_out=js["summary"]["after_filtering"]["total_reads"]
            r1_in,r1_kept = parse_tg(f"results/trimgalore/{s}/R1_trimming_report.txt")
            r2_in,r2_kept = parse_tg(f"results/trimgalore/{s}/R2_trimming_report.txt")
            mean=median=mode=pct=None
            with open(f"results/bbmerge/{s}/merge_stats_{s}.txt") as ih:
                for ln in ih:
                    if ln.startswith('#PercentOfPairs'):
                        pct=float(ln.strip().split()[-1])
                    if ln.startswith('#Mean'):
                        mean=float(ln.strip().split()[-1])
                    if ln.startswith('#Median'):
                        median=float(ln.strip().split()[-1])
                    if ln.startswith('#Mode'):
                        mode=float(ln.strip().split()[-1])
            rows.append([s,fp_in,fp_out,r1_in,r1_kept,r2_in,r2_kept,pct,mean,median,mode])
        with open(output[0], 'w', newline='') as out:
            w=csv.writer(out, delimiter='\t')
            w.writerows(rows)



# --- Notebook execution (relative abundance) ---
rule run_relative_abundance_nb:
    conda: "envs/analysis.yaml"
    input:
        table="results/qiime/export/feature-table-frequency.tsv",
        tax="results/qiime/export/taxonomy/taxonomy.tsv"
    params:
        nb=lambda w: config["analysis"]["notebook"]
    output:
        nb="results/analysis/relative_abundance.executed.ipynb",
        html="results/analysis/relative_abundance.html"
    shell:
        r"""
        mkdir -p results/analysis
        # Execute the notebook in-place into results/analysis
        jupyter nbconvert \
          --to notebook --execute "{params.nb}" \
          --output "relative_abundance.executed" \
          --output-dir results/analysis \
          --ExecutePreprocessor.timeout=1200

        # Convert to HTML for easy viewing
        jupyter nbconvert \
          --to html results/analysis/relative_abundance.executed.ipynb \
          --output relative_abundance.html \
          --output-dir results/analysis
        """
# --- Notebook execution (relative abundance) ---
rule run_relative_abundance_nb:
    conda: "envs/analysis.yaml"
    input:
        table="results/qiime/export/feature-table-frequency.tsv",
        tax="results/qiime/export/taxonomy/taxonomy.tsv"
    params:
        nb=lambda w: config["analysis"]["notebook"]  # e.g., notebooks/relative_abundance.ipynb
    output:
        nb="results/analysis/relative_abundance.executed.ipynb"
    shell:
        r"""
        mkdir -p results/analysis
        jupyter nbconvert \
          --to notebook --execute "{params.nb}" \
          --output "relative_abundance.executed" \
          --output-dir results/analysis \
          --ExecutePreprocessor.timeout=1200 \
          --ExecutePreprocessor.kernel_name=python3
        """
