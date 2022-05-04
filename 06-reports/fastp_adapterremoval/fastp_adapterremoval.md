Can fastp replace AR2 or leeHOM for adapter removal and overlapping
merged reads?
================
Alex Huebner
May 04, 2022

-   <a href="#results" id="toc-results">Results</a>
    -   <a href="#efficiency-to-remove-adapter-sequences"
        id="toc-efficiency-to-remove-adapter-sequences">Efficiency to remove
        adapter sequences</a>
    -   <a href="#merging-overlapping-sequencing-pairs-into-a-single-molecule"
        id="toc-merging-overlapping-sequencing-pairs-into-a-single-molecule">Merging
        overlapping sequencing pairs into a single molecule</a>
-   <a href="#conclusion" id="toc-conclusion">Conclusion</a>
-   <a href="#references" id="toc-references">References</a>

The trimming of adapter sequences is the first step of processing
Illumina sequencing data that is required to obtain clean data for
subsequent analyses. Due to the short average DNA molecule length,
paired-end sequencing of ancient DNA leads to an either complete or
partial overlap of the mate sequences. Therefore, specialised adapter
removal tools, such as AdapterRemoval2 (AR2; Schubert, Lindgreen, and
Orlando (2016)) or leeHOM (Renaud, Stenzel, and Kelso 2014) were
developed to allow the merging of these read pairs into a single
sequence and became the de-facto standard in the ancient DNA field.

Both these programs operate on FastQ files. Due to the recent change in
the sequence chemistry of the Illumina sequencers, which introduced
poly-G strings at the end of short DNA sequences, an additional
processing step is required to remove these artefacts from the actual
biological sequences. The tool that is commonly applied for this step is
fastp (Chen et al. 2018) and is a very fast, general purpose FastQ file
manipulation tool. In a recent version, the authors added the option to
trim adapter sequences and merge overlapping read pairs, too. Therefore,
the question arose whether it would be necessary to still run the
specialised adapter removal programs or whether fastp can replace these
and save an additional processing step.

In the following, I will compare the three programs, AdapterRemoval2,
leeHOM, and fastp, with each other on a simulated set of sequencing data
that varied regarding the amount of ancient DNA damage and the average
DNA molecule length.

# Results

To evaluate the performance of the different tools, I made use of a
previously simulated dataset (for details see
<http://github.com/alexhnr/aDNA-DenovoAssembly>) that consisted out of
50 million DNA sequences that were generated from 30 microbial species
commonly found in dental calculus as either commensals or contaminants.
The samples differed in average DNA molecule length (short, medium,
long) and amount of observed ancient DNA damage (none, light, heavy).
Therefore, I tested each of the three programs across these nine samples
using the default settings, enabling or disabling the collapse of
overlapping read pairs, respectively. I summarised the performance of
the adapter removal using FastQC (Andrews et al. 2010) and summarised
all reports into summary tables using multiQC (Ewels et al. 2016).

## Efficiency to remove adapter sequences

First and for most, the importance of the adapter removal step is to
remove the technical adapter sequences from the sequencing data because
they might interfere with the analysis of the biological signal in these
data. Therefore, I first evaluated whether the softwares could remove
the adapter sequences efficiently and no overrepresented adapter
sequences were observed in the data. For this, I used FastQC (Andrews et
al. 2010) and summarised it using multiQC (Ewels et al. 2016).

The first metric to evaluate is whether the respective sequencing data
of the samples still contains sequences that are likely Illumina adapter
sequences. FastQC iterates over the sequencing data searching for such
sequence motifs and reports a sample either as pass, warn, or fail.
Since FastQC is running its analyses per sequencing data file, I further
summarised this metric across all output files of a sample.

I compared the results when either without (**Figure 1**) or with
(**Figure 2**) the merging of overlapping read pairs. Without merging
the read pairs, both AdapterRemoval2 and fastp returned one output file
per input file (in total two files), while leeHOM additionally created a
single-end file in case only one mate of the read pairs passed the
required quality. For all nine samples, both AdapterRemoval2 and fastp
managed to pass the adapter content test of FastQC, i.e. no adapter
sequence motifs were observed. For leeHOM, there were no adapter
sequencing motifs observed in the singleton data but the forward mate of
the read pairs failed the test when the underlying DNA molecule length
distribution was short. Such a DNA molecule length distribution is
characteristic for poorly preserved ancient DNA samples and is
particularly worrying, when planning to use leeHOM to remove adapter
sequences prior to assembly.

![](fastp_adapterremoval_files/figure-gfm/fastqc_adaptercontent_pe-1.png)<!-- -->

**Figure 1**: **Overview of the results of the analysis for the absence
of Illumina adapter sequence motifs in the sequencing data after the
adapter removal without merging overlapping reads.** The fraction of
FastQ files was dependent on the number of files reported per program:
both AdapterRemoval2 and fastp exported each one file per read mate,
while leeHOM exported additionally singleton for which one of the read
mates did not fulfil the quality criteria.

When enabling the merging of overlapping reads, the overlapping and
subsequently merged sequencing data exported by all three programs
passed the adapter content test of FastQC. However, for the read pairs
that could not be overlapped we observed differences between the
programs. While fastp had no adapter sequence motifs in the remaining
paired data independent of the underlying DNA molecule length
distribution, we observed issues for both AdapterRemoval2 and leeHOM,
when we had either medium or short read length distribution.

![](fastp_adapterremoval_files/figure-gfm/fastqc_adaptercontent_se-1.png)<!-- -->

**Figure 2**: **Overview of the results of the analysis for the absence
of Illumina adapter sequence motifs in the sequencing data after the
adapter removal with merging overlapping reads.** The fraction of FastQ
files was dependent on the number of files reported per program.

Neither program retained any other overrepresented sequences and all
three programs raised either warnings or even failed the per-base
sequence content test of FastQC, when the underlying DNA molecule
distribution was short, suggesting that this is rather a problem with
the underlying FastQC test rather than something that these programs
influenced.

In summary, fastp performed the most consistent across all data sets
regarding the removal of adapter sequences. While AdapterRemoval2 only
seemed to have issues regarding removing adapters in non-overlapping
sequence pairs when the DNA molecule length distribution was shorter,
leeHOM seemed to suffer in general with removing adapter sequences in
sequences that were not merged into a single one.

## Merging overlapping sequencing pairs into a single molecule

Ancient DNA molecules can be characterised by their on-average short
length that often is shorter than the nominal sequencing length of the
sequencing kit. For such cases, the sequences obtained from paired-end
sequencing overlap in the middle and can be merged into a single, longer
DNA molecule sequence. This is commonly performed for all type of
analyses with the exception of *de novo* assembly, which either requires
or strongly benefits from paired-end sequencing data. Therefore, I
additionally will evaluate the performance of merging overlapping paired
reads.

The behaviour of merging of overlapping read pairs can usually be
influenced by specifying the required overlap of the two paired reads.
By default, fastp requires an overlap of 30 bp, AdapterRemoval2 of 11
bp, and for leeHOM this parameter is not accessible via the commandline.
The effect of this can be observed when comparing the fraction of merged
read pairs per program (**Figure 3**). AdapterRemoval2 had the highest
fraction, followed by leeHOM and fastp, however, the differences were
relatively small, particularly for the samples simulated from a medium
or short DNA molecule length distribution.

![](fastp_adapterremoval_files/figure-gfm/fraction_merged-1.png)<!-- -->

**Figure 3**: **Overview of the fraction of merged overlapping read
pairs per program.**

However, while the number of merged read pairs was highest for
AdapterRemoval2, the average read length of the merged molecules was
much shorter than leeHOM and fastp, which returned almost identical
sequencing length (**Figure 4**). This suggests that enabling the
quality trimming options `--trimns --trimqualities` that we enable by
default substantially trimmed the read pairs while the other two
programs do not perform such extensive quality trimming by default.

![](fastp_adapterremoval_files/figure-gfm/avg_length-1.png)<!-- -->

**Figure 4**: **Average read length of merged overlapping read pairs per
program.**

# Conclusion

To evaluate whether the program fastp is a suitable alternative for
trimming Illumina adapter sequences from sequencing data, I compared it
to two standard adapter removal tools, AdapterRemoval2 and leeHOM. I
selected a set of simulated sequencing data with different
characteristics regarding the DNA molecule length and the amount of
ancient DNA damage and ran all programs with parameters commonly used in
the field of ancient DNA. Fastp was able to remove adapter sequences
efficiently and it was the only program that could remove all adapter
sequence motifs and could thus pass the adapter content test of FastQC.
Furthermore, it is able to efficiently merge overlapping read pairs into
a single sequence. Overall, it could therefore be used to replace the
dedicated adapter trimming tools, AdapterRemoval2 and leeHOM, because it
provides an efficient alternative to them.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Andrews2010" class="csl-entry">

Andrews, Simon et al. 2010. “FastQC: A Quality Control Tool for High
Throughput Sequence Data.” Babraham Bioinformatics, Babraham Institute,
Cambridge, United Kingdom.

</div>

<div id="ref-Chen2018" class="csl-entry">

Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. 2018. “Fastp: An
Ultra-Fast All-in-One FASTQ Preprocessor.” *Bioinformatics* 34 (17):
i884–90. <https://doi.org/10.1093/bioinformatics/bty560>.

</div>

<div id="ref-Ewels2016" class="csl-entry">

Ewels, Philip, Måns Magnusson, Sverker Lundin, and Max Käller. 2016.
“MultiQC: Summarize Analysis Results for Multiple Tools and Samples in a
Single Report.” *Bioinformatics* 32 (19): 3047–48.
<https://doi.org/10.1093/bioinformatics/btw354>.

</div>

<div id="ref-RenaudLeeHOM2014" class="csl-entry">

Renaud, Gabriel, Udo Stenzel, and Janet Kelso. 2014. “<span
class="nocase">leeHom</span>: Adaptor Trimming and Merging for Illumina
Sequencing Reads.” *Nucleic Acids Research* 42 (18): e141.
<https://doi.org/10.1093/nar/gku699>.

</div>

<div id="ref-Schubert2016" class="csl-entry">

Schubert, Mikkel, Stinus Lindgreen, and Ludovic Orlando. 2016.
“AdapterRemoval V2: Rapid Adapter Trimming, Identification, and Read
Merging.” *BMC Research Notes* 9 (1): 88.
<https://doi.org/10.1186/s13104-016-1900-2>.

</div>

</div>
