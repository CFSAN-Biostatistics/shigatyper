# ShigaTyper

ShigaTyper is a quick and easy tool designed to determine Shigella serotype using Illumina paired end reads with low computation requirement. It has been tested on an Ubuntu 16.04.3 LTS guest addition on VMware Player version 14.1.1 in a Windows 7 and Virtual Box version 5.2.4 in a Windows 10 operating system.

## Usage

    usage: shigatyper [-h] [-n SAMPLE_NAME] [--verbose] [--version] read1 read2

    positional arguments:
    read1
    read2

    optional arguments:
    -h, --help            show this help message and exit
    -n SAMPLE_NAME, --name SAMPLE_NAME
    --verbose, -v
    --version             show program's version number and exit

## Example

    $ shigatyper CFSAN029786_S10_L001_R1_001.fastq.gz CFSAN029786_S10_L001_R2_001.fastq.gz
    sample  prediction      ipaB
    CFSAN029786     Shigella dysenteriae serotype 3 +