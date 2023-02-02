################################################################################
# Project: Tiny projects
# Subproject: Sequence annotation with Bakta
#
# In a recent tweet, the author of Prokka, Torsten Seeman, recommended to check
# the annotation pipeline Bakta (https://github.com/oschwengers/bakta) because
# it is still actively maintained compared to Prokka. Prokka is the
# quasi-standard for annotating metagenome-assembled genomes but comes with some
# issues, next to the lack of database updates, e.g. the very slow step to
# create GenBank files.
#
# In the following, I will compare the annotations of Bakta to the ones obtained
# by Prokka in order to see whether Bakta would be indeed a valid candidate for
# replacing Prokka in the future.
#
# Alex Huebner, 12/10/22
################################################################################

from glob import glob
import os

import numpy as np
import pandas as pd
import pyfastx

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
MAGS, = glob_wildcards("03-data/emn001_mags/{mag}.fasta.gz")
################################################################################

rule all:
    input:
        "05-results/002_Bakta_EMN001_annotations.tsv.gz"

#### Bakta #####################################################################

rule bakta:
    output:
        "04-analysis/002_mag_annotation/bakta/{mag}.tsv"
    message: "Annotate with Bakta: {wildcards.mag}"
    conda: "002_bakta.yaml"
    resources:
        mem = 32,
        cores = 8 
    params:
        fasta = "03-data/emn001_mags/{mag}.fasta.gz",
        dir = "04-analysis/002_mag_annotation/bakta",
        db = "/mnt/archgen/users/huebner/refdbs/bakta/db"
    threads: 8
    shell:
        """
        bakta -p {wildcards.mag} \
            --db {params.db} \
            -o {params.dir} \
            --keep-contig-headers \
            --threads {threads} \
            {params.fasta}
        """

rule summarise_bakta:
    input:
        expand("04-analysis/002_mag_annotation/bakta/{mag}.tsv", mag= MAGS)
    output:
        "05-results/002_Bakta_EMN001_annotations.tsv.gz"
    message: "Summarise the results of Bakta into a single table"
    run:
        reports = [pd.read_csv(fn, sep="\t", skiprows=2) \
                     .assign(binID=os.path.basename(fn).replace(".tsv", ""))
                   for fn in input]

        pd.concat(reports) \
            .rename({'#Sequence Id': "contigID"}, axis=1) \
            .iloc[:, [9] + list(range(9))] \
            .to_csv(output[0], sep="\t", index=False, compression="gzip")

################################################################################

#### Prokka ####################################################################

rule summarise_prokka:
    output:
        "05-results/002_Prokka_EMN001_annotations.tsv.gz"
    message: "Summarise the results of Prokka into a single table"
    params:
        dir = "/mnt/archgen/microbiome_paleobiotech/EMN001_Paleofuran/04-analysis/automatic_MAG_refinement/aDNA_samples_human/EMN001-megahit/bins"
    run:
        def read_gff_to_tbl(fn):
            tbl_lines = []
            with open(fn, "rt") as gfffile:
                for line in gfffile:
                    if line.rstrip() == "##FASTA":
                        break
                    else:
                        if line[:2] != "##":
                            tbl_lines.append(tuple(line.rstrip().split("\t")))
            return tbl_lines

        gffs = []
        for fn in glob(f"{params.dir}/*.gff"):
            gff = pd.DataFrame(read_gff_to_tbl(fn),
                        columns=['contig_name', 'source', 'Type', 'Start', 'Stop',
                                'score', 'Strand', 'frame', 'DbXrefs']) \
                .drop(['source', 'score', 'frame'], axis=1) \
                .assign(binID=os.path.basename(fn).replace(".gff", "")) \
                .query('Type != "gene"')
            orig_contig_names = {prokka[0]: orig[0]
                                 for prokka, orig in zip(pyfastx.Fasta(fn.replace("gff", "fna"), build_index=False),
                                                         pyfastx.Fasta(fn.replace("gff", "fasta.gz"), build_index=False))}
            gff['contig_name'] = [orig_contig_names[c] for c in gff['contig_name'].values]
            gffs.append(gff)

        pd.concat(gffs) \
            .iloc[:, [6] + list(range(6))] \
            .sort_values(['binID', 'contig_name']) \
            .to_csv(output[0], sep="\t", index=False, compression="gzip")

################################################################################


