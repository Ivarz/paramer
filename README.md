# Paramer
`paramer` is a k-mer based sequence search tool designed to search for parasite sequences within shotgun metagenomic datasets.
By leveraging Bloom filter and nthash hashing algorithm, `paramer` can quickly identify and classify parasitic DNA from sequencing samples.
The tool is proof of concept and needs additional testing to be considered production ready.

## Dependencies
* zlib
* cmake
* ntHash 2.2.0

## Installation

```
git clone --recurse-submodules https://github.com/Ivarz/paramer
cd paramer
git -C third_party/ntHash checkout 2.2.0
```

## Compilation

To compile:
```
make
```

To perform unit tests:
```
make test
```


## Usage
```
Usage:
  paramer <COMMAND>

  where command is:
      mask           mask fasta file with kraken2 file and/or kmers found in another fasta file
      bloom-build    build Bloom's filter
      bloom-search   search in Bloom's filter
      extend         extend sequence in 3' and 5' directions with kmers found in Bloom's filter
      stats-bloom    gather metrics from bloom filter
      stats-fasta    gather metrics from fasta files
```

Tool can utilize output from [kraken2](https://github.com/DerrickWood/kraken2) to perform reference sequence masking.

Full parasite analysis workflow consists of following steps:

1. acquire reference sequences of interest (e.g reference genomes of _Ascaris lumbricoides_)
```
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ascaris_lumbricoides/PRJEB4950/ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.fa.gz
```
2. perform reference classification with [kraken2](https://github.com/DerrickWood/kraken2). 
For gut microbiome sequencing datasets I recommend using host specific kraken2 catalogue (if available).
Assuming [kraken2](https://github.com/DerrickWood/kraken2) is available in `$PATH`:

```
kraken2 \
        --threads 24 \
        --output k2.out.txt \
        --report k2.rep.txt \
        --db /path/to/kraken2/database \
        ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.fa.gz
```

3. perform reference masking from the obtained [kraken2](https://github.com/DerrickWood/kraken2) results:
```
paramer mask \
    -f ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.fa.gz \
    -k k2.out.txt \
    -o ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.masked.fa.gz
```

4. build bloom filter from the reference
```
paramer bloom-build \
    -k 31 \
    -w 31 \
    -s 10G \
    -r ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.masked.fa.gz \
    -o ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.masked.blm
```

5. search your sequencing data against Bloom filter.
Assuming paired-end sequencing sample names `sample_1.fq.gz` and `sample_2.fq.gz`:
```
paramer bloom-search \
        -b ascaris_lumbricoides.PRJEB4950.WBPS19.genomic.masked.blm \
        -i sample_1.fq.gz \
        -I sample_2.fq.gz \
        -o sample_output_1.fq.gz \
        -O sample_output_2.fq.gz
```

Will output files containing sequencing reads that have matches in the Bloom filter 
(in example: `sample_output_1.fq.gz` and `sample_output_2.fq.gz`)
