# danbing-tk
A toolkit for variable number tandem repeats (VNTRs) analysis, which enables:
1. building repeat-pan genome graphs (RPGG) for a set of "genotypable" VNTRs given haplotype-resolved assemblies and matched short-read seqeuncing (SRS) data (`danbing-tk build`),
2. genotyping each VNTR as a set of (*k*-mer, count) given SRS data (`danbing-tk align`), and
3. predicting VNTR length from the genotype (`danbing-tk predict`) given locus-specific sampling bias (LSB).


## Installation
This will take up to an hour on a typical work station, depending on your internet connection.

### Download Releases
Each release comes with the lastest version of genotypable VNTRs, RPGG, and precomputed LSB.

|                   | File                              | Input of           | Output of        |
|-------------------|-----------------------------------|--------------------|------------------|
| Genotypable VNTRs | tr.good.bed                       | danbing-tk build   |                  |
| RPGG              | pan.(tr\|lntr\|rntr\|graph).kmers | danbing-tk align   | danbing-tk build |
| precomputed LSB   | step1.csv                         | danbing-tk predict |                  |

- Genotypable VNTRs are also available under`data/tr.good.bed` if RPGG and LSB are not needed.
- The 397 non-repetitive control regions are available under `data/ctrl.bed`.

### Install Dependencies
danbing-tk requires several external packages. It is recommended to install all requirements using conda as follows:

```bash
conda install -c conda-forge -c bioconda -c intel \
    snakemake=5.11.2 samtools=1.10 bedtools=2.29.2 minimap2=2.17 \
    scikit-learn=0.23.1 matplotlib=3.3.0
```

If the requirements are in conflict with existing packages, create a new environment specifically for danbing-tk with:
```
conda create -n $MY_ENVIRONMENT -c conda-forge -c bioconda -c intel \
    snakemake=5.11.2 samtools=1.10 bedtools=2.29.2 minimap2=2.17 \
    scikit-learn=0.23.1 matplotlib=3.3.0
```

### Building on Linux
```shell
git clone https://github.com/ChaissonLab/danbing-tk
cd danbing-tk && make
```

### Test Environment
To check if everything is configured properly:
1. Go to `/$PREFIX/danbing-tk/test/`
2. Replace `$PREFIX` in `goodPanGenomeGraph.json` and `input/bam.fofn` with the path to danbing-tk
3. Run `snakemake -p -s ../pipeline/GoodPanGenomeGraph.snakefile -j 4 --rerun-incomplete --output-wait 3`

Tested on V0.0. 
## Usage

### danbing-tk align
Decompress the RPGG `RPGG.tar.gz` and link `*.kmers` (required by the `-qs` option) to your working directory with `ln -s`.

An example usage to genotype SRS sample using the RPGG:

```shell
samtools fasta -@2 -n $SRS.bam |
bam2pe -fai /dev/stdin |
danbing-tk -gc 50 -k 21 -qs pan -fai /dev/stdin -o $OUT_PREFIX \
           -p 24 -cth 45 -rth 0.5
```

`danbing-tk align` takes ~42 cpu hours to genotype a 30x SRS sample. This will generate a `*.tr.kmers` output with format specified in [File Format](#file_format).

**Important note:** If outputs of `danbing-tk align` are intended to be used directly for downstream analyses e.g. association tests, please check the [distribution of LSB](#distribution-of-lsb) section below before running.

### danbing-tk build
Required inputs: haplotype-resolved assemblies, matched SRS data, reference genome (major chromosomes only without minor contigs) and bed file of tandem repeat regions (optional)

Copy `/$PREFIX/danbing-tk/pipeline/goodPanGenomeGraph.snakefile` to your working directory and edit accordingly. Run the snakemake pipline with:

```bash
snakemake -p -s /$PREFIX/danbing-tk/pipeline/GoodPanGenomeGraph.snakefile -j 40\
    --cluster "{params.copts} -c {resources.cores} --mem={resources.mem}G -k" \
    --rerun-incomplete --restart-times 1 --output-wait 30
```

Submitting jobs to cluster is preferred as `danbing-tk build` is compute-intensive, ~1200 cpu hours for the original dataset. Otherwise, remove `--cluster` and its parameters to run jobs locally.

### danbing-tk predict
Locus-specific sampling biases (LSB) at VNTR regions are critical for normalizing the sum of *k*-mer counts to VNTR length. We provided precomputed LSB at the VNTR regions for fast comparison, however this assumes the LSB of the dataset of interest is close enough to the dataset in the original paper. Please ensure this assumption is valid by running a joint PCA on the LSB of non-repetitive regions with the original dataset, provided in `ctrl.cov`. If this assumption failed, leave-one-out analysis (next section) on the dataset of interest is necessary to make accurate predictions. The following usage is for when the assumption holds.

Link precomputed statistics of the original dataset to your working directory.

`ln -s /$PREFIX/danbing-tk/dat/step1.csv /$WORKING_DIR/analysis/.`

Run length prediction with:

```bash
/$PREFIX/danbing-tk/script/kmc2length.py --genome GENOME.TXT --nloci NLOCI \
    --skip1 --LOOconf LOOCONF --sampleConf SAMPLECONF
```

Predictions and accuracies are written to `pred_len.txt` and `rel_err.txt`

## Analysis
### Distribution of LSB
To ensure technical variations are consistent between samples and do not introduce bias to VNTR genotypes, or *k*-mer counts, in a subset of samples. It is necessary to ensure the distribution of LSB shows up in one cluster. An example analysis is shown in [`script/LSB_analysis.ipynb`](https://github.com/ChaissonLab/danbing-tk/blob/master/script/LSB_analysis.ipynb) and `script/getCovByLocus.397.sh`.

### Leave-one-out analysis
To evaluate the quality of custom RPGG on matching SRS dataset, copy `/$PREFIX/danbing-tk/pipeline/leaveOneOut.snakefile` to your working directory and edit accordingly.
Run the snakemake pipleine with:

```bash
snakemake -p -s /$PREFIX/pipeline/LeaveOneOut.snakefile -j 40 --cluster \
    "{params.copts} -c {resources.cores} --mem={resources.mem}G -k" \
    --rerun-incomplete --restart-times 1 --output-wait 30
```

Submitting jobs to cluster is preferred as `danbing-tk build` is compute-intensive; otherwise, remove `--cluster` and its parameters to run jobs locally.


## File Format
- \*.graph.kmers
```
>locus i
kmer0	out_edges0
kmer1	out_edges1
...
>locus i+1
...
```
`out_edges` denotes the presence of T/G/C/A as the next nucleotide encoded with 4 bits.

- \*.(tr|ntr).kmers
```
>locus i
kmer0	kmer_count0
kmer1	kmer_count1
...
>locus i+1
...
```
