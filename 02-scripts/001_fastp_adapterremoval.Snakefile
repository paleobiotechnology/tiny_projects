################################################################################
# Project: Tiny projects
# Subproject: Adapter removal with fastp
#
# Fastp is a versatile program that includes a number of different quality
# filtering steps, such as discarding low-quality reads, trimming
# poly-nucleotide strings at read ends etc. One of these steps also includes the
# identification and removal of adapter sequences.
#
# So far, the ancient DNA community has commonly applied the programs
# AdapterRemoval2 [Schubert2013] and leeHOM [Renaud2015] because these programs
# can merge overlapping reads. However, this option is not unique to these two
# tools but can also be performed by fastp.
#
# Here, we want to test the capabilities of adapter removal by the three
# programs fastp, AdapterRemoval2, and leehom with and without overlapping
# merged reads using simulated double-stranded library preparation, paired-end
# sequencing data that were generated for the aDNA metagenome assembly project.
#
# Alex Huebner, 14/03/22
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### Sequencing data ###########################################################
SAMPLES, = glob_wildcards("/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_1.fq.gz")
MERGING = ['collapse', 'nocollape']
TOOLS = ['adapterremoval', 'leeHOM', 'fastp']
################################################################################

rule all:
    input:
        general = "05-results/001_fastp_adapterremoval-general_stats.tsv",
        fastqc = "05-results/001_fastp_adapterremoval-fastqc.tsv"

#### AdapterRemoval2 ###########################################################

rule adapterremoval_summary:
    input: 
        "04-analysis/001_fastp_adapterremoval/adapterremoval_multiqc.html"
    output:
        touch("04-analysis/001_fastp_adapterremoval/adapterremoval.done")
    params:
        dir = "04-analysis/001_fastp_adapterremoval/adapterremoval"
    shell:
        """
        find {params.dir} -mindepth 1 -maxdepth 1 -type f \( ! -name "*.log" -and ! -name "*.done" \) -exec rm {{}} +
        """

rule adapterremoval:
    output:
        touch("04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}.done")
    message: "Remove adapter sequences using AdapterRemoval2: {wildcards.sample} with option {wildcards.merging}"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 16,
        cores = 8
    params:
        pe1 = "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_1.fq.gz",
        pe2 = "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_2.fq.gz",
        prefix = "04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}",
        forward_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGACTATCTCGTATGCCGTCTTCTGCTTG",
        reverse_adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        collapse = lambda wildcards: "--minadapteroverlap 1 --collapse" if wildcards.merging == "collapse" else ""
    log: "04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}.adapterremoval.log"
    threads: 8
    shell:
        """
        AdapterRemoval --basename {params.prefix} \
                           --file1 {params.pe1} \
                           --file2 {params.pe2} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --threads {threads} \
                           --qualitybase 33 \
                           --adapter1 {params.forward_adapter} \
                           --adapter2 {params.reverse_adapter} \
                           {params.collapse} \
                           --settings {log}
        """

rule fastqc_adapterremoval:
    input:
        "04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}.done"
    output:
        touch("04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}.fastqc") 
    message: "Run FastQC on the FastQ files of Adpaterremoval: {wildcards.sample} with option {wildcards.merging}"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        prefix = "04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}",
    shell:
        """
        mkdir -p {params.prefix}/
        fastqc -o {params.prefix} \
               -t {threads} \
               {params.prefix}.{{pair1.truncated,pair2.truncated,singleton.truncated,collapsed,collapsed.truncated}}
        """

rule multiqc_adapterremoval:
    input:
        expand("04-analysis/001_fastp_adapterremoval/adapterremoval/{sample}-{merging}.fastqc", sample=SAMPLES, merging=MERGING)
    output:
        "04-analysis/001_fastp_adapterremoval/adapterremoval_multiqc.html"
    message: "Summarise adapterremoval using multiQC"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        prefix = "04-analysis/001_fastp_adapterremoval/adapterremoval",
    shell:
        """
        multiqc -d {params.prefix} -f -o multiqc && \
        multiqc -d {params.prefix} -f -n {output}
        """

################################################################################

### leeHOM #####################################################################


rule leeHOM_summary:
    input: 
        "04-analysis/001_fastp_adapterremoval/leeHOM_multiqc.html"
    output:
        touch("04-analysis/001_fastp_adapterremoval/leeHOM.done")
    params:
        dir = "04-analysis/001_fastp_adapterremoval/leeHOM"
    shell:
        """
        find {params.dir} -mindepth 1 -name "*.fq.gz" -exec rm {{}} +
        """

rule leeHOM:
    output:
        touch("04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}.done")
    message: "Remove adapter sequences using leeHom: {wildcards.sample} with option {wildcards.merging}"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 16,
        cores = 8
    params:
        pe1 = "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_1.fq.gz",
        pe2 = "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_2.fq.gz",
        prefix = "04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}",
        forward_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGACTATCTCGTATGCCGTCTTCTGCTTG",
        reverse_adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        collapse = lambda wildcards: "--ancientdna " if wildcards.merging == "collapse" else ""
    log: "04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}.leeHOM.log"
    threads: 8
    shell:
        """
        leeHom -f {params.forward_adapter} \
               -s {params.reverse_adapter} \
               -t {threads} \
               -fq1 {params.pe1} \
               -fq2 {params.pe2} \
               -fqo {params.prefix} \
               --log {log} \
               {params.collapse}
        """

rule fastqc_leeHOM:
    input:
        "04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}.done"
    output:
        touch("04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}.fastqc") 
    message: "Run FastQC on the FastQ files of leeHOM: {wildcards.sample} with option {wildcards.merging}"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        prefix = "04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}",
    shell:
        """
        mkdir -p {params.prefix}/
        fastqc -o {params.prefix} \
               -t {threads} \
               {params.prefix}{{_r1,_r2,}}.fq.gz
        """

rule multiqc_leeHOM:
    input:
        expand("04-analysis/001_fastp_adapterremoval/leeHOM/{sample}-{merging}.fastqc", sample=SAMPLES, merging=MERGING)
    output:
        "04-analysis/001_fastp_adapterremoval/leeHOM_multiqc.html"
    message: "Summarise leeHOM using multiQC"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        prefix = "04-analysis/001_fastp_adapterremoval/leeHOM",
    shell:
        """
        multiqc -d {params.prefix} -f -o multiqc && \
        multiqc -d {params.prefix} -f -n {output}
        """


################################################################################

#### FastP #####################################################################

rule fastp_summary:
    input:
        "04-analysis/001_fastp_adapterremoval/fastp_multiqc.html"
    output:
        touch("04-analysis/001_fastp_adapterremoval/fastp.done")
    params:
        dir = "04-analysis/001_fastp_adapterremoval/leeHOM"
    shell:
        """
        find {params.dir} -mindepth 1 -name "*.fastq.gz" -exec rm {{}} +
        """


rule fastp:
    output:
        touch("04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}.done")
    message: "Remove adapter sequences using FastP: {wildcards.sample} with option {wildcards.merging}"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 16,
        cores = 8
    params:
        pe1 = "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_1.fq.gz",
        pe2 = "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/simulated_samples/50M/50M_{sample}_0_2.fq.gz",
        prefix = "04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}",
        forward_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGACTATCTCGTATGCCGTCTTCTGCTTG",
        reverse_adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        collapse = lambda wildcards: f"--merge --merged_out 04-analysis/001_fastp_adapterremoval/fastp/{wildcards.sample}-{wildcards.merging}_collapsed.fastq.gz" if wildcards.merging == "collapse" else ""
    log: "04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}.fastp.json"
    threads: 8
    shell:
        """
        fastp --in1 {params.pe1} \
              --in2 {params.pe2} \
              --out1 {params.prefix}_1.fastq.gz \
              --out2 {params.prefix}_2.fastq.gz \
              {params.collapse} \
              --compression 4 \
              --adapter_sequence {params.forward_adapter} \
              --adapter_sequence_r2 {params.reverse_adapter} \
              --trim_poly_g \
              --json {log} \
              --html /dev/null
        """

rule fastqc_fastp:
    input:
        "04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}.done"
    output:
        touch("04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}.fastqc") 
    message: "Run FastQC on the FastQ files of Fastp: {wildcards.sample} with option {wildcards.merging}"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        prefix = "04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}",
    shell:
        """
        mkdir -p {params.prefix}/
        fastqc -o {params.prefix} \
               -t {threads} \
               {params.prefix}{{_1,_2,_collapsed}}.fastq.gz
        """

rule multiqc_fastp:
    input:
        expand("04-analysis/001_fastp_adapterremoval/fastp/{sample}-{merging}.fastqc", sample=SAMPLES, merging=MERGING)
    output:
        "04-analysis/001_fastp_adapterremoval/fastp_multiqc.html"
    message: "Summarise fastp using multiQC"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        prefix = "04-analysis/001_fastp_adapterremoval/fastp",
    shell:
        """
        multiqc -d {params.prefix} -f -o multiqc && \
        multiqc -d {params.prefix} -f -n {output}
        """

################################################################################

#### Summarise the MultiQC reports #############################################

rule multiqc_summary:
    input:
        ar2 = "04-analysis/001_fastp_adapterremoval/adapterremoval_multiqc.html",
        leeHOM = "04-analysis/001_fastp_adapterremoval/leeHOM_multiqc.html",
        fastp = "04-analysis/001_fastp_adapterremoval/fastp_multiqc.html"
    output:
        general = "05-results/001_fastp_adapterremoval-general_stats.tsv",
        fastqc = "05-results/001_fastp_adapterremoval-fastqc.tsv"
    message: "Summarise the output tables generated by MultiQC"
    conda: "001_fastp_adapterremoval.yaml"
    resources:
        mem = 8
    params:
        dir = "04-analysis/001_fastp_adapterremoval"
    script:
        "rscripts/001_fastp_adapterremoval-multiqc_summary.R"

################################################################################
