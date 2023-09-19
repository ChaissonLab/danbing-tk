# danbing-tk
A toolkit for variable number tandem repeats (VNTRs) analysis, which enables:
1. building repeat-pan genome graphs (RPGG) given haplotype-resolved assemblies for genome-wide profiling or simply VNTR alleles for targeted genotyping (referred to as `danbing-tk build`),
2. genotyping each VNTR as a set of (*k*-mer, count) given short-read sequencing (SRS) data (referred to as `danbing-tk align`), and
3. estimating VNTR or motif dosage from the genotype with bias correction (referred to as `danbing-tk predict`).

Our [initial manuscript](https://www.nature.com/articles/s41467-021-24378-0) illustrates the key concept of this tool.
The [latest update](https://genome.cshlp.org/content/33/4/511) details improvements on QC and bias correction, and extended applications on eQTL discoveries with motif compositions.


## Download Releases
The latest release [v1.3.1-manuscript](https://github.com/ChaissonLab/danbing-tk/releases/tag/v1.3.1-manuscript) comes with the lastest version of VNTR set, RPGG, and QC statistics.

|                   | File                                                                   | Input of             | Output of                                 |
|-------------------|------------------------------------------------------------------------|----------------------|-------------------------------------------|
| VNTR set          | `tr.good.bed`                                                          | `danbing-tk build`   |                                           |
| RPGG              | `pan.tr.kmers`, `pan.kmerDBi.umap`, `pan.kmerDBi.vv`, `pan.graph.umap` | `danbing-tk align`   | `danbing-tk build` or `vntr2kmers_thread` |

- Release v1.3.1-manuscript: provided QC statistics and all resources associated with the latest manuscript.
- Release v1.3: Updated RPGG built from 35 HGSVC genomes.
- Release v1.0: VNTR summary statistics and eGene discoveries are also included. Example analyses such as differential length/motif analysis, eQTL mapping, VNTR locus QC, sample QC are also included.


## Building on Linux
```shell
git clone --recursive https://github.com/ChaissonLab/danbing-tk
cd danbing-tk && make -j 5
```

## danbing-tk align
Decompress the RPGG `RPGG.tar.gz` in your working directory.

An example usage to genotype SRS sample using the RPGG:

```shell
samtools fasta -@2 -n $SRS.bam |
/$PREFIX/danbing-tk/bin/danbing-tk -gc 85 3 -ae -kf 4 1 -cth 45 -o $OUT_PREF -k 21 -qs pan -fa /dev/stdin -p $THREADS | gzip >$OUT_PREF.aln.gz
```

`danbing-tk align` takes ~12 cpu hours to genotype a 30x SRS sample. This will generate `$OUT_PREF.tr.kmers` and `$OUT_PREF.aln.gz` output with format specified in [File Format](#file-format).

**Important note:** If outputs of `danbing-tk align` are intended to be compared across individuals e.g. association studies, please check the bias_correction [notebook](https://github.com/ChaissonLab/eMotif_manuscript_analysis_scripts/tree/main/bias_correction) before running.


## danbing-tk build
### Install Dependencies
For users intend to use `danbing-tk align` or the [Scenario 1](#scenario-1-building-an-rpgg-for-a-single-tr-locus-given-vntr-alleles) of `danbing-tk build`, this step is not required.

The `danbing-tk build` pipeline and `danbing-tk predict` require several external packages. It is recommended to install all requirements using conda as follows:

```bash
conda create -n $MY_ENVIRONMENT -c conda-forge -c bioconda \
    python=3.11.4 snakemake=7.30.1 minimap2=2.26 samtools=1.17 bedtools=2.31.0 statsmodels=0.14.0 matplotlib=3.7.2
conda activate $MY_ENVIRONMENT
```

To check if everything is configured properly (tested on v1.3.2):
1. Go to `/$PREFIX/danbing-tk/test/`
2. Replace `$PREFIX` in `goodPanGenomeGraph.json` and `input/genome.bam.tsv` with the path to danbing-tk
3. Run `snakemake -p -s ../pipeline/GoodPanGenomeGraph.snakefile -j 4 --forceall --output-wait 3`


### Running danbing-tk build
#### Scenario 1: building an RPGG for a single TR locus given VNTR alleles
- Required inputs:
	- VNTR alleles for each haplotype (one FASTA per haplotype)

- Run `vntr2kmers_thread` with something like this:
	- `vntr2kmers_thread -g -k 21 -fs 700 -ntr 700 -on $NAME -fa $NUM_HAPLOTYPES $LIST_OF_FASTAS`
	- **Note**: at least 500 bp flanks are required for accurate mapping of pair-end reads, 700 bp was specified in the above example.

- Index the graph as follows to use `danbing-tk align` later:
	- `/$PREFIX/danbing-tk/bin/ktools serialize $NAME`


#### Scenario 2: building an RPGG for a VNTR set given assemblies

- Required inputs: 
	- haplotype-resolved assemblies (FASTA)
	- matched SRS data (BAM; optional)
	- GRCh38 (FASTA; major chromosomes only without minor contigs)
	- tandem repeat regions (BED; `tr.good.bed` from the [release page](https://github.com/ChaissonLab/danbing-tk/releases/tag/v1.3.1-manuscript) or user-defined)

- Copy `/$PREFIX/danbing-tk/pipeline/goodPanGenomeGraph.json` to your working directory and edit accordingly. 

- A config file `/$INPUT_DIR/genome.bam.tsv` with two columns, one for genome name and one for bam file path, is required if SRS data is available for graph pruning, e.g.

    ```HG00514 /panfs/qcb-panasas/tsungyul/HG00514/HG00514.IL.srt.bam```
    
    Otherwise, set `pruning` in `goodPanGenomeGraph.json` to `False` and use a single column input for `genome.bam.tsv`.

- Run the snakemake pipline with:
```bash
snakemake -p -s /$PREFIX/danbing-tk/pipeline/GoodPanGenomeGraph.snakefile -j 40\
    --cluster "{params.copts} -c {resources.cores} --mem={resources.mem}G -k" \
    --rerun-incomplete --restart-times 1 --output-wait 30
```

Submitting jobs to cluster is preferred as `danbing-tk build` is compute-intensive, ~1200 cpu hours for the original dataset. Otherwise, remove `--cluster` and its parameters to run jobs locally.


## danbing-tk predict
Locus- and sample-specific biases are critical for normalizing the sum of *k*-mer counts to VNTR dosage (as a proxy for predicted length) and normalizing the average of *k*-mer counts to motif dosage. The bias for each locus in each sample is computed from the deviation of local read depth from the global mean given a set of *invariant k-mers*. Examples of this analysis can be found [here](https://github.com/ChaissonLab/eMotif_manuscript_analysis_scripts/tree/main/bias_correction) 

Automated bias correction has been added since [v1.3.2](https://github.com/ChaissonLab/danbing-tk/releases/tag/v1.3.2). Invariant *k*-mer metadata `ikmer.meta` (human readable version `ikmer.meta.txt`) and example `trkmers.meta.txt` can be found in [Assets](https://github.com/ChaissonLab/danbing-tk/releases/tag/v1.3.2).

Example usage:
```
/$PREFIX/danbing-tk/bin/danbing-tk-pred trkmers.meta.txt ikmer.meta corrected.gt.tsv bias.tsv
```

Caveat: Estimated *k*-mer dosage could be inaccurate if the bias term is too close to zero.

## Miscellaneous
### Leave-one-out analysis
To evaluate the quality of custom RPGG on matching SRS dataset, copy `/$PREFIX/danbing-tk/pipeline/leaveOneOut.snakefile` to your working directory and edit accordingly.
Run the snakemake pipleine with:

```bash
snakemake -p -s /$PREFIX/pipeline/LeaveOneOut.snakefile -j 40 --cluster \
    "{params.copts} -c {resources.cores} --mem={resources.mem}G -k" \
    --rerun-incomplete --restart-times 1 --output-wait 30
```

Submitting jobs to cluster is preferred as this analysis is compute-intensive; otherwise, remove `--cluster` and its parameters to run jobs locally.


## File Format
### \*.graph.kmers
```
>locus i
kmer0	out_edges0
kmer1	out_edges1
...
>locus i+1
...
```
`out_edges` denotes the presence of T/G/C/A as the next nucleotide encoded with 4 bits.

### \*.(tr|ntr).kmers
```
>locus i
kmer0	kmer_count0
kmer1	kmer_count1
...
>locus i+1
...
```
The second field is optional.

**Important Note**: the output of `danbing-tk align` do not contain locus info and the first field for minimal disk usage. The table can be reconstructed using the `danbing_aln_output.tr_kmers.metadata.txt.gz` from `metadata.tar.gz` on [Zenodo](https://sandbox.zenodo.org/record/1169833#.ZAo3FNLMKEJ)

### Alignment output (`-a` option)
- Synopsis
	```
	<src> <dest> <read_name> <read_seq/0> <read_seq/1> <ops/0> <annot/0> <ops/1> <annot/1>
	```
- `src`: source locus of a read pair (for simulation only)
- `dest`: aligned locus for the read pair
- `ops`: nucleotide-level operations to align the read to the graph. `size = read_len + #del`
	- `=`: a match 
	- `X[A|C|G|T]`: a mismatch; letter in the graph is shown
	- `D[A|C|G|T]`: a deletion in the read; letter in the graph is shown
	- `I`: an insertion in the read
	- `H`: a nucleotide in the homopolymer run; only shown when length > ksize
	- `*`: unalinged nucleotide
- `annot`: kmer-level VNTR annotations after applying `ops` to the read. `size = read_len - ksize + 1 + #del - #ins`
	- `=`: a match in the repeat
	- `.`: a match in the flank
	- `*`: unaligned kmer
