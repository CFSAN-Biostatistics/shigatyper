# ShigaTyper

ShigaTyper is a quick and easy tool designed to determine Shigella serotype using Illumina or Oxford Nanopore reads
with low computation requirement.

## Usage

```{bash}
usage: shigatyper.py [-h] [--R1 FASTA] [--R2 FASTA] [--SE FASTA] [--ont] [-n SAMPLE_NAME] [--verbose] [--version]

ShigaTyper v. 2.0.0, 2022

options:
  -h, --help            show this help message and exit
  --R1 FASTA            Input FASTQ is R1 of paired-end reads
  --R2 FASTA            Input FASTQ is R2 of paired-end reads
  --SE FASTA            Input FASTQ is contains single-end reads
  --ont                 The input FASTQ file contains ONT reads
  -n SAMPLE_NAME, --name SAMPLE_NAME
  --verbose, -v
  --version             show program's version number and exit
```

## Example

```{bash}
# Paired-end reads
shigatyper.py --R1 SRX5006488_R1.fastq.gz --R2 SRX5006488_R2.fastq.gz
sample  prediction      ipaB
SRX5006488      Shigella boydii serotype 12     +

# Single-end reads
shigatyper.py --SE SRX5006488.fastq.gz
sample  prediction      ipaB
SRX5006488-se   Shigella boydii serotype 12     +

# Oxford Nanopore reads
shigatyper.py --SE SRX7050861.fastq.gz --ont
sample  prediction      ipaB
SRX7050861-ont  Shigella dysenteriae serotype 3 +
```
