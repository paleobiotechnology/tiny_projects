library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# Load files
## General stats
general_stats <- map(c("adapterremoval", "leeHOM", "fastp"), function(p)
                     fread(str_c(snakemake@params[["dir"]], p,
                                 "multiqc_data/multiqc_general_stats.txt", sep = "/")))
names(general_stats) <- c("adapterremoval", "leeHOM", "fastp")

## FastQC
fastqc <- map(c("adapterremoval", "leeHOM", "fastp"), function(p)
              fread(str_c(snakemake@params[["dir"]], p,
                          "multiqc_data/multiqc_fastqc.txt", sep = "/")))
names(fastqc) <- c("adapterremoval", "leeHOM", "fastp")

# Process and summarise data
## General stats
gs_summary <- list()
gs_summary[["adapterremoval"]] <- general_stats[["adapterremoval"]] %>%
  mutate(program = "adapterremoval",
         settings = str_split_fixed(Sample, " \\| ", n=5)[,4]) %>%
  separate(settings, c("rlength", "damage", "collapsed"), sep = "_|-") %>%
  mutate(collapsed = collapsed == "collapse",
         seqtype = if_else(str_detect(Sample, "collapsed") | str_detect(Sample, "singleton"),
                           "SE", "PE")) %>%
  filter(!str_detect(Sample, "pair2")) %>%
  select(program:collapsed, seqtype,
         avg_length = `FastQC_mqc-generalstats-fastqc-avg_sequence_length`,
         total_sequence = `FastQC_mqc-generalstats-fastqc-total_sequences`) %>%
  group_by(program, rlength, damage, collapsed, seqtype) %>%
  summarise(avg_length = mean(avg_length),
            total_sequence = sum(total_sequence)) %>%
  ungroup()
gs_summary[["leeHOM"]] <- general_stats[["leeHOM"]] %>%
  filter(!str_detect(Sample, "leeHOM$")) %>%
  mutate(program = "leeHOM",
         settings = str_split_fixed(Sample, " \\| ", n=5)[,4]) %>%
  separate(settings, c("rlength", "damage", "collapsed"), sep = "_|-") %>%
  mutate(collapsed = collapsed == "collapse",
         seqtype = if_else(str_detect(Sample, "collapsed") | str_detect(Sample, "singleton"),
                           "SE", "PE")) %>%
  filter(!str_detect(Sample, "_r2")) %>%
  select(program:collapsed, seqtype,
         avg_length = `FastQC_mqc-generalstats-fastqc-avg_sequence_length`,
         total_sequence = `FastQC_mqc-generalstats-fastqc-total_sequences`) %>%
  group_by(program, rlength, damage, collapsed, seqtype) %>%
  summarise(avg_length = mean(avg_length),
            total_sequence = sum(total_sequence)) %>%
  ungroup()
gs_summary[["fastp"]] <- general_stats[["fastp"]] %>%
  filter(!str_detect(Sample, "-collapse$")) %>%
  mutate(program = "fastp",
         settings = str_split_fixed(Sample, " \\| ", n=5)[,4]) %>%
  separate(settings, c("rlength", "damage", "collapsed"), sep = "_|-") %>%
  mutate(collapsed = collapsed == "collapse",
         seqtype = if_else(str_detect(Sample, "collapsed") | str_detect(Sample, "singleton"),
                           "SE", "PE")) %>%
  filter(!str_detect(Sample, "_2")) %>%
  select(program:collapsed, seqtype,
         avg_length = `FastQC_mqc-generalstats-fastqc-avg_sequence_length`,
         total_sequence = `FastQC_mqc-generalstats-fastqc-total_sequences`) %>%
  group_by(program, rlength, damage, collapsed, seqtype) %>%
  summarise(avg_length = mean(avg_length, na.rm = T),
            total_sequence = sum(total_sequence, na.rm = T)) %>%
  ungroup()
gs_summary <- bind_rows(gs_summary)

## FastQC
fqc_summary <- list()
fqc_summary[["adapterremoval"]] <- fastqc[["adapterremoval"]] %>%
  mutate(program = "adapterremoval",
         settings = str_split_fixed(Sample, " \\| ", n=5)[,4]) %>%
  separate(settings, c("rlength", "damage", "collapsed"), sep = "_|-") %>%
  mutate(collapsed = collapsed == "collapse",
         seqtype = if_else(str_detect(Sample, "collapsed") | str_detect(Sample, "singleton"),
                           "SE", "PE")) %>%
  filter(!str_detect(Sample, "singleton")) %>%
  select(program:collapsed, seqtype,
         per_base_sequence_content, overrepresented_sequences, adapter_content) %>%
  mutate(across(per_base_sequence_content:adapter_content, factor, levels = c("pass", "warn", "fail")))
fqc_summary[["leeHOM"]] <- fastqc[["leeHOM"]] %>%
  mutate(program = "leeHOM",
         settings = str_split_fixed(Sample, " \\| ", n=5)[,4]) %>%
  separate(settings, c("rlength", "damage", "collapsed"), sep = "_|-") %>%
  mutate(collapsed = collapsed == "collapse",
         seqtype = if_else(str_detect(Sample, "collapsed") | str_detect(Sample, "singleton"),
                           "SE", "PE")) %>%
  select(program:collapsed, seqtype,
         per_base_sequence_content, overrepresented_sequences, adapter_content) %>%
  mutate(across(per_base_sequence_content:adapter_content, factor, levels = c("pass", "warn", "fail")))
fqc_summary[["fastp"]] <- fastqc[["fastp"]] %>%
  mutate(program = "fastp",
         settings = str_split_fixed(Sample, " \\| ", n=5)[,4]) %>%
  separate(settings, c("rlength", "damage", "collapsed"), sep = "_|-") %>%
  mutate(collapsed = collapsed == "collapse",
         seqtype = if_else(str_detect(Sample, "collapsed") | str_detect(Sample, "singleton"),
                           "SE", "PE")) %>%
  select(program:collapsed, seqtype,
         per_base_sequence_content, overrepresented_sequences, adapter_content) %>%
  mutate(across(per_base_sequence_content:adapter_content, factor, levels = c("pass", "warn", "fail")))
fqc_summary <- bind_rows(fqc_summary)

# Write files
fwrite(gs_summary, sep = "\t", file = snakemake@output[["general"]])
fwrite(fqc_summary, sep = "\t", file = snakemake@output[["fastqc"]])
