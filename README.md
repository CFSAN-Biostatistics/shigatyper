# ShigaTyper

ShigaTyper is a quick and easy tool designed to determine Shigella serotype using Illumina (single or paired-end) or Oxford Nanopore reads
with low computation requirement.

## Installation

Shigatyper is available from [Bioconda](https://bioconda.github.io/recipes/shigatyper/README.html) and can be installed using
the following command.

```{bash}
conda create -n shigatpyer -c conda-forge -c bioconda shigatyper
```

## Running ShigaTyper

ShigaTyper supports compressed FASTQs as inputs. These FASTQs can be single-end or paired-end Illumina reads, or reads from
Oxford Nanopore.

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

## Example Runs

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

## ShigaTyper Outputs

After your run is complete, two tab-delimited TSV files (`<PREFIX>.tsv` and `<PREFIX>-hits.tsv`) are created with the results. By default 
the output files, uses the base name of the input FASTQ file. You can change this by using the `--name` parameter.

### Example `<PREFIX>.tsv`

This file contains the final serotype predicted by ShitaTyper. It looks like the following:

```{bash}
# With A predicted Serotype
sample	prediction	ipaB	notes
SRX7050861-ont	Shigella dysenteriae serotype 3	+	this strain is ipaB+, suggesting that it retains the virulent invasion plasmid.

# Note Shigella or EIEC
sample	prediction	ipaB	notes
ERR3772599	Not Shigella or EIEC	-	No read was mapped to the reference sequence database.
```

The `<PREFIX>.tsv` will have the following four collumns.

| Column Name | Description                                 |
|-------------|---------------------------------------------|
| sample      | The name of the input sample                |
| prediction  | The serotype predicted by ShigaTyper        |
| ipaB        | The precence of ipaB (`+`) or absence (`-`) |
| notes       | Any notes associated with result            |

### Example `<PREFIX>-hits.tsv`

The `<PREFIX>-hits.tsv` file will contain statistics about each individual gene hit. If there are no hits, this file
will not be produced (e.g. non-Shigella or EIEC inputs).

Here's an example of how it will look:

```{bash}
	Hit	Number of reads	Length Covered	reference length	% covered	Number of variants	% accuracy
0	ipaH_c	331	780	780	100.0	10	98.7
1	ipaB	59	1743	1743	100.0	44	97.5
2	Sd3_wzx	18	1515	1515	100.0	7	99.5
3	Sd3_wzy	20	1104	1104	100.0	3	99.7
```

| Column Name        | Description                                            |
|--------------------|--------------------------------------------------------|
| index              | Index number in the array                              |
| Hit                | Name of the gene                                       |
| Number of reads    | Number of reads mapped to the Hit                      |
| Length Covered     | Length of reference gene aligned to                    |
| reference length   | Length of the reference gene                           |
| % covered          | Percent of the reference gene aligned to               |
| Number of variants | Number of varaints in the alignment                    |
| % accuracy         | Percent of identical matches across the reference gene |

## Citations

If you make use of this tool, please cite the following:

* **[ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper)**  
This tool, for serotyping Shigella.
Wu Y, Lau HK, Lee T, Lau DK, Payne J [In Silico Serotyping Based on Whole-Genome Sequencing Improves the Accuracy of Shigella Identification.](https://doi.org/10.1128/AEM.00165-19) *Applied and Environmental Microbiology*, 85(7). (2019)  

* **[BCFtools](https://github.com/samtools/bcftools)**  
Utilities for variant calling and manipulating VCFs and BCFs.  
Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H [Twelve years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) *GigaScience* Volume 10, Issue 2 (2021)  

* **[Minimap2](https://github.com/lh3/minimap2)**  
A versatile pairwise aligner for genomic and spliced nucleotide sequences  
Li H [Minimap2: pairwise alignment for nucleotide sequences.](https://doi.org/10.1093/bioinformatics/bty191) *Bioinformatics* 34:3094-3100 (2018)

* **[Samtools](https://github.com/samtools/samtools)**  
Tools for manipulating next-generation sequencing data  
Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R [The Sequence Alignment/Map format and SAMtools](http://dx.doi.org/10.1093/bioinformatics/btp352). *Bioinformatics* 25, 2078–2079 (2009)  
