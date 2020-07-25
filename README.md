## danbing-tk
A toolkit for variable number tandem repeats (VNTRs) analysis, which enables:
1. building repeat-pan genome graphs for a set of "genotypable" VNTRs given haplotype-resolved assemblies (**danbing-tk build**), 
2. genotyping each VNTR as a set of (*k*-mer, count) given short-read sequencing (SRS) data (**danbing-tk align**), and
3. predicting VNTR length from the genotype (**danbing-tk predict**).


### Compiling source code

`g++ -std=c++11 -O3 -lpthead danbing-tk.cpp`

`g++ -std=c++11 -O3 bam2pe.cpp`

- repeat-pan genome graph databse

    `./kmers/pan.*.kmers` are used by **danbing-tk align**

- sample usage

```shell
samtools fasta -@2 -n srs.bam |
awk '{if (substr($1,1,1) == ">") {
        if (substr($1,length($1)-1,1) == "/") { print substr($1, 1, length($1)-2) } else { print $1 } }
      else { print $1 }
     }' |
bam2pe -k 21 -fai /dev/stdin |
danbing-tk -gc 50 -k 21 -qs pan -fai /dev/stdin -o output_prefix -p 24 -cth 45 -rth 0.5
```

- Output format

```
>locus i
kmer0    kmer_count0
kmer1    kmer_count1
...
>locus i+1
...
```
