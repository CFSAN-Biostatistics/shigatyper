#!/usr/bin/env python3
# coding: utf-8
"""
usage: shigatyper.py [-h] [--R1 FASTA] [--R2 FASTA] [--SE FASTA] [--ont] [-n SAMPLE_NAME] [--verbose] [--version]

A WGS-based genoserotyping pipeline for Shigella spp.

options:
  -h, --help            show this help message and exit
  --R1 FASTA            Input FASTQ is R1 of paired-end reads
  --R2 FASTA            Input FASTQ is R2 of paired-end reads
  --SE FASTA            Input FASTQ is contains single-end reads
  --ont                 The input FASTQ file contains ONT reads
  -n SAMPLE_NAME, --name SAMPLE_NAME
  --verbose, -v
  --version             show program's version number and exit
"""

from . import __version__
version = __version__ # in case I didn't find every usage

usage = f"""ShigaTyper v. {__version__}, 2022

A WGS-based genoserotyping pipeline for Shigella spp.

Yun Wu, Henry K Lau, Teresa Lee, David K Lau, Justin Payne

    The bacteria Shigella spp., consisting of 4 species and >50
serotypes, cause shigellosis, a foodborne disease of significant
morbidity, mortality, and economic loss worldwide. Classical Shigella
identification based on selective media and serology is tedious,
time-consuming, expensive, and not always accurate. Molecular diagnostic
assay does not distinguish Shigella at species level or from
enteroinvasive Escherichia coli (EIEC). We inspected the whole genome
sequencing (WGS) data from 219 Shigella isolates and observed low
concordance rate between conventional designation and molecular
serotyping, 86.8% and 78.9% at species and serotype level, respectively.
Serotype determinants for 6 additional serotypes were identified.
Examination of differentiation gene markers commonly perceived as
characteristic hallmarks in Shigella showed high variability among
different serotypes. Using this information, we developed ShigaTyper, an
automated workflow that utilizes limited computational resources to
accurately and rapidly determine 58 Shigella serotypes using Illumina
paired end WGS reads. Shigella serotype determinants and species-specific
diagnostic markers were first identified through read alignment to an
in-house curated reference sequence database. Relying on sequence hits
that passed a threshold level of coverage and accuracy, serotype can be
unambiguously predicted within 1 min for an average sized WGS sample of
~500 MB. Validation with WGS data from 112 isolates show an accuracy of
98.2%. This pipeline is the first step towards building a comprehensive
WGS-based analysis pipeline of Shigella spp. in a field laboratory
setting, where speed is essential and resources need to be more
cost-effectively dedicated.

"""
import os
import csv
import datetime
import logging
import itertools
import pandas as pd
import shutil
import sys
import tempfile

from subprocess import check_output, CalledProcessError
from os.path import join as j, dirname

rlog = logging.getLogger('ShigaTyper') #root log

def pairwise(i, n=2):
    args = [iter(i)] * n
    return itertools.zip_longest(*args)

def readable(tdelta):
    one_milli = datetime.timedelta(milliseconds=1)
    one_second = datetime.timedelta(seconds=1)
    one_minute = datetime.timedelta(minutes=1)
    one_hour = datetime.timedelta(hours=1)
    if tdelta <= one_milli:
        return f"{tdelta.microseconds} microseconds"
    elif one_milli < tdelta <= one_second:
        return f"{tdelta.microseconds / 1000 : 02} milliseconds"
    elif one_second < tdelta <= one_minute:
        return f"{tdelta.seconds : 02} seconds"
    elif one_minute < tdelta <= one_hour:
        return f"{tdelta.seconds // 60 } minutes, {tdelta.seconds} seconds"
    else:
        return f"{tdelta.seconds // 3600} hours, {tdelta.seconds // 60 : 01} seconds"

def run(reads, tempdir, sample_name='', threshold=50, rlog=rlog, ont=False, outdir='./', *args, **kwargs):
    # ## 3. Map Filtered reads to Reference sequence database (ShigellaRef5)
    rlog.critical(f"ShigaTyper v{__version__}")

    rlog = rlog.getChild('run')
    
    timetrack = []
    def sub(cmd, log=rlog):
        log = log.getChild(cmd.split()[0]) #log with name of subprocess command
        result = str(check_output(cmd, shell=True, encoding='UTF-8')).strip()
        log.info(cmd)
        for row in result.splitlines():
            log.debug("{}".format(row.replace(' ', '\t')))
        return result

    read1 = os.path.abspath(reads[0])
    read2 = None

    if len(reads) == 2:
        read2 = os.path.abspath(reads[1])

    log = rlog.getChild('index')
    start = datetime.datetime.now()
    # give an error message if the reference sequence database is not there
    ShigellaRef = os.path.abspath(j(dirname(__file__), 'resources', 'ShigellaRef5.fasta'))
    if not os.path.isfile(ShigellaRef):
        log.error(f"Error: reference sequence not found at {ShigellaRef}")
        exit()
    # indexing reference sequence database when necessary
    mmi_index = j(tempdir, "ShigellaRef5.mmi")
    if not os.path.isfile(mmi_index):
        log.warning(f"building Reference sequence index in {mmi_index}")
        sub(f'minimap2 -d {mmi_index} {ShigellaRef}', log)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    log = rlog.getChild('map')
    start = datetime.datetime.now()
    # map the fastq.gz files to reference sequence database
    outputbam = j(tempdir, "_Shigella5.bam")

    if read2:
        # paired-end reads
        sub(f'minimap2 -ax sr {mmi_index} {read1} {read2} | samtools view -F 0x04 -b | samtools sort -o {outputbam}', log)
    else:
        # single-end reads
        if ont:
            # ont reads
            sub(f'minimap2 -ax map-ont {mmi_index} {read1} | samtools view -F 0x04 -b | samtools sort -o {outputbam}', log)
        else:
            sub(f'minimap2 -ax sr {mmi_index} {read1} | samtools view -F 0x04 -b | samtools sort -o {outputbam}', log)

    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ### Checkpoint 1 
    # If there is no read mapped to the reference sequence database, discontinue the analysis.
    log = rlog.getChild('check1')
    start = datetime.datetime.now()
    checkpoint = 0
    ipaB = 0
    prediction = ""
    NR = sub(f'samtools view {outputbam} | wc -l', log)
    if int(NR[0]) == 0:
        checkpoint = 1
        log.error("Checkpoint 1 failed.")
        prediction ="Not Shigella or EIEC"
    else: log.info("Checkpoint 1......passed.")
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ## 3.1. Examine sequence hits identified by Minimap2
    log = rlog.getChild('hits')
    start = datetime.datetime.now()
    #check what sequences were hits and how many reads were mapped to each of the hits
    if checkpoint == 0:
        Hits = sub(f'samtools view {outputbam} | cut -f3 | uniq -c', log)
        hits = []
        Nreads = []
        # 527 ipaB
        # 299 Sd3_wzx
        # 150 Sd3_wzy
        for nread, hit in pairwise(Hits.split()):
            hits.append(hit)
            Nreads.append(int(nread))
        Maphits = pd.DataFrame({'Hit': hits, 'Number of reads': Nreads})
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ### Checkpoint 2
    # If the strain is not Shigella or EIEC, discontinue the analysis.
    log = rlog.getChild('check2')
    start = datetime.datetime.now()
    if checkpoint == 0:
        if "ipaH_c" not in hits: 
            checkpoint = 2
            if "Sb13_wzx" in hits:
                prediction = "Shigella boydii serotype 13"
                checkpoint = 13
            else:
                log.error("Checkpoint 2 failed!")
                prediction = "Not Shigella or EIEC"
        else: log.info("Checkpoint 2.....passed.")
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ## 3.2. Determine breadth of coverage and accuracy of the mapped hits
    log = rlog.getChild('accuracy')
    start = datetime.datetime.now()
    if checkpoint == 0:
        # find reference sequence length for calculation of % coverage
        Gene_length = list()
        RefDic = dict()
        Reflines = sub(f'samtools view -H {outputbam} | grep "SN"', log)
        # @SQ     SN:Sb12_wzx     LN:1335
        # @SQ     SN:Sb12_wzy     LN:1140
        # @SQ     SN:Sb13_wzx     LN:1320
        for _, sn, ln in pairwise(Reflines.split(), 3):
            sn = sn.replace('SN:', '')
            ln = ln.replace('LN:', '')
            RefDic[sn] = int(ln)
        for target, length in RefDic.items():
            if target in hits:
                Gene_length.append(length)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())


    log = rlog.getChild('pileup')
    start = datetime.datetime.now()
    if checkpoint == 0:
        # index the bam file for mpileup
        sub(f'samtools index {outputbam}', log)
        # index reference sequence for mpileup if the index is not already there
        outputmpileup = j(tempdir, "Shigella5.mpileup")
        if ont:
            sub(f'samtools mpileup -C10 -B -f {ShigellaRef} {outputbam} -o {outputmpileup}', log)
        else:
            sub(f'samtools mpileup -C50 -q 20 -Q 20 -f {ShigellaRef} {outputbam} -o {outputmpileup}', log)
        # can I pipe it so that I won't have a mpileup file left in the disk?
        lapse = datetime.datetime.now() - start
        log.info(f"Complete in {readable(lapse)}.")
        timetrack.append(lapse.total_seconds())
    else: log.info("skipped.")

    log = rlog.getChild('hits')
    start = datetime.datetime.now()
    if checkpoint == 0:
        covSummary = sub(f"cat {outputmpileup} | awk '$4 > 0' | cut -f1 | uniq -c", log)
        hits = []; bpCovered = []
        # 1654 ipaB
        # 1495 Sd3_wzx
        for cover, hit in pairwise(covSummary.split()):
            hits.append(hit)
            bpCovered.append(int(cover))
        hitsCovered = pd.DataFrame({'Hit': hits, 'Length Covered': bpCovered})
        # merge 2 tables on "Hit"
        List1 = pd.merge(left=Maphits, right = hitsCovered, left_on="Hit", right_on = "Hit", how = 'left')
        List1 = List1.fillna(0)
        List1["reference length"] = Gene_length
        # calculate percent reference sequence covered
        List1["% covered"] = round(100* List1["Length Covered"].div(List1["reference length"]), 1)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    log = rlog.getChild('variants')
    start = datetime.datetime.now()
    if checkpoint == 0:
        variants = None
        if ont:
            variants = sub(f'samtools mpileup -C10 -B -f {ShigellaRef} -g {outputbam} | bcftools call -m | cat | grep -v "^#" | grep PL | cut -f1 | uniq -c', log)
        else:
            variants = sub(f'samtools mpileup -C50 -q 20 -Q 20 -f {ShigellaRef} -g {outputbam} | bcftools call -m | cat | grep -v "^#" | grep PL | cut -f1 | uniq -c', log)
        hits = []; Nvar = []
        # 7 ipaH_c
        # 17 ipaB
        for nvar, hit in pairwise(variants.split()):
            hits.append(hit)
            Nvar.append(int(nvar))
        Variants = pd.DataFrame({'Hit':hits, 'Number of variants':Nvar})

    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    log = rlog.getChild('re-axis')
    start = datetime.datetime.now()
    if checkpoint == 0:
        # merge the two tables and add % accuracy
        List2 = pd.merge(left=List1, right = Variants, left_on="Hit", right_on = "Hit", how = 'left')
        List2 = List2.fillna(0)
        List2["accurate seq"] = List2["Length Covered"] - List2["Number of variants"]
        List2["% accuracy"] = round(100*List2["accurate seq"].div(List2["Length Covered"]), 1)
        List2 = List2.drop('accurate seq', axis =1)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    log = rlog.getChild('filter')
    start = datetime.datetime.now()
    if checkpoint == 0:
        # filter based on threshold for % coverage and 80% for % accuracy
        FList = List2[(List2["% covered"] >= threshold) & (List2["% accuracy"] >= 80)]
        log.info("Analysis completed.")
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ### Checkpoint 3 
    # If the coverage for the sample is too low, or the strain is lacY+ (potentially EIEC), discontinue the analysis.
    log = rlog.getChild('check3')
    start = datetime.datetime.now()
    if checkpoint == 0:
        Hits = FList.Hit.tolist()
        if "ipaH_c" in Hits:
            Hits.remove("ipaH_c")
        else:
            checkpoint = 31
            if "Sb13_wzx" in hits:
                prediction = "Shigella boydii serotype 13"
                checkpoint = 13
            else:
                log.error("Checkpoint 3 failed!")
                prediction = "Not Shigella or EIEC"
            
        # ipaB may have lower coverage because of its locating on a plasmid though.
        if "ipaB" in Hits:
            ipaB += 1
            Hits.remove("ipaB")
        
        if "EclacY" in Hits:
            exception = ['Sb9_wzx', 'Sb15_wzx']
            if any(gene in Hits for gene in exception) is True:
                checkpoint = 0
            else:
                checkpoint = 32
                prediction = "EIEC"
    if checkpoint == 0: log.info("Checkpoint 3.....passed.")
    elif checkpoint == 3: log.error("Checkpoint 3 failed!")
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ### Checkpoint 4 
    # If there are no hits at all left, or are multiple hits of O-antigen serotype determinant, discontinue the analysis.
    log = rlog.getChild('check4')
    start = datetime.datetime.now()
    if checkpoint == 0:
        if len(Hits) == 0:
            checkpoint = 41
            prediction = "No prediction (no wzx)"
        else:
            def wzxindex(Hits):
                ind = []
                for i, j in enumerate(Hits):
                    if j.find("wzx") > 0: ind.append(i)
                return ind

            if len(wzxindex(Hits)) > 1:
                checkpoint = 42
                prediction = "No prediction (multiple wzx)"
    else: log.info("Skipped.")
                
    if checkpoint == 0: log.info("Checkpoint 4 ..... passed")

    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # ## 4. Shigella serotype prediction
    log = rlog.getChild('predict')
    start = datetime.datetime.now()
    if checkpoint == 0:
        if "cadA" in Hits:
            Hits.remove("cadA")
            if "Ss_methylase" in Hits:
                if len(wzxindex(Hits)) == 0:
                    if ipaB == 0: 
                        prediction = "Shigella sonnei form II"
                    else: 
                        prediction = "Shigella sonnei (low levels of form I)"
                elif "Ss_wzx" in Hits:
                    prediction = "Shigella sonnei, form I"
                else:
                    prediction = "EIEC"; checkpoint = 32
            elif "Sd1_wzx" in Hits:
                if "Sd1_rfp" in Hits: 
                    prediction = "Shigella dysenteriae serotype 1"
                else: 
                    prediction = "Shigella dysenteriae serotype 1, rfp- (phenotypically negative)"
            elif "Sd8_wzx" in Hits:
                prediction = "Shigella dysenteriae serotype 8"
            elif "Sb11_wzx" in Hits:
                prediction = "Shigella boydii serotype 11"
            else:
                prediction = "EIEC"; checkpoint = 32
        else:
        ## the rationale is the S. boydii 6 has poor quality at the IS insertion junction (bp 252:253) so the ratio of 
        ## mpileup depth to samtools depth at the junction is way lower than the ratio for the overall gene
            if "Sb6_wzx" in Hits:
                if "wbaM" in Hits:
                    a = sub(f'samtools depth -r wbaM:252-253 {outputbam} | cut -f3')
                    a = a.split("\n")
                    # From here the script is failing. Because there is a \n in variable a(also in b,c,d down from here).
                    # Splitting on \n makes variable a(and b,c,d) an iterable list of items.
                    # When making this adjustment for variable a,b,c and d, the pipeline returns the correct serotype.
                    a = [int(i) for i in a]; average_a = sum(a)/len(a)
                    if a == 0:
                        prediction = "Shigella boydii serotype 6 or 10"
                    else:
                        b = sub(f"cat {outputmpileup} | awk '$1==\"wbaM\" && $2 >251 && $2 <254' | cut -f4")
                        b = b.split("\n")
                        b = [int(i) for i in b]; average_b = sum(b)/len(b)
                        junction_ratio = average_b/average_a
                        c = sub(f'samtools depth -r wbaM {outputbam} | cut -f3')
                        c = c.split("\n")
                        c = [int(i) for i in c]; average_c = sum(c)/len(c)
                        d = sub(f'cat {outputmpileup} | awk \'$1=="wbaM"\' | cut -f4')
                        d = d.split("\n")
                        d = [int(i) for i in d]; average_d = sum(d)/len(d)
                        overall_ratio = average_d/average_c        
                        if junction_ratio/overall_ratio > 0.5:
                            prediction = "Shigella boydii serotype 10"
                        else:
                            prediction = "Shigella boydii serotype 6"
                else: prediction = "Shigella boydii serotype 6 or 10" 
            else:
                if len(wzxindex(Hits)) == 0:
                    prediction = "No prediction (no wzx)"
                    checkpoint = 41
                else: 
                    wzx = Hits[wzxindex(Hits)[0]]
                    if wzx == "Sb1_wzx":
                        if "heparinase" in Hits:
                            prediction = "Shigella boydii serotype 20"
                        else: 
                            prediction = "Shigella boydii serotype 1"
                    elif wzx == "SbProv_wzx":
                        prediction = "Shigella boydii Provisional serotype E1621-54"
                    elif wzx == "SdProv_wzx":
                        prediction = "Shigella dysenteriae Provisional serotype 96-265"
                    elif wzx == "SdProvE_wzx":
                        prediction = "Shigella dysenteriae Provisional serotype E670-74"
                    elif wzx[1] == "b":
                        prediction = "Shigella boydii serotype " + wzx[2:wzx.find("_")]
                    elif wzx[1] == "d":
                        prediction = "Shigella dysenteriae serotype " + wzx[2:wzx.find("_")]
                    elif wzx == "Sf6_wzx":
                        prediction = "Shigella flexneri serotype 6"
                    else:
                        try: Hits.remove("Sf_wzy")
                        except: Hits = Hits
                        if len(Hits) == 0:
                            prediction = "Shigella flexneri serotype Y"
                        else:
                            SfDic = {"Shigella flexneri Yv": ["Xv"], "Shigella flexneri serotype 1a": ["gtrI"], 
                                "Shigella flexneri serotype 1b": ["gtrI", "Oac1b"], "Shigella flexneri serotype 2a": 
                                ["gtrII"], 'Shigella flexneri 2av': ['gtrII', 'Xv'], 
                                "Shigella flexneri serotype 2b": ["gtrII", "gtrX"], "Shigella flexneri serotype 3a":
                                ["gtrX","Oac"], "Shigella flexneri serotype 3b": ["Oac"], "Shigella flexneri serotype 4a": 
                                ["gtrIV"], "Shigella flexneri serotype 4av": ["gtrIV", "Xv"], "Shigella flexneri serotype 4b":
                                ["gtrIV", "Oac"], 'Shigella flexneri 4bv': ['gtrIV', 'Oac', 'Xv'],
                                "Shigella flexneri serotype 5a": (["gtrV", "Oac"], ['gtrV']),
                                    "Shigella flexneri serotype 5b": (["gtrV", "gtrX", "Oac"], ['gtrV', 'gtrX']),
                                    "Shigella flexneri serotype X": ["gtrX"], "Shigella flexneri serotype Xv (4c)":
                                ["gtrX", "Xv"], "Shigella flexneri serotype 1c (7a)": ['gtrI', 'gtrIC'], 
                                "Shigella flexneri serotype 7b": ['gtrI', "gtrIC", "Oac1b"]}
                            predict = 0
                            for Serotype, Targets in SfDic.items():
                                if Targets == Hits:
                                    prediction = Serotype; predict += 1
                                elif Hits in Targets:
                                    prediction = Serotype; predict += 1
                            if predict == 0:
                                prediction = "Shigella flexneri, novel serotype"
        lapse = datetime.datetime.now() - start
        log.info(f"Complete in {readable(lapse)}.")
        timetrack.append(lapse.total_seconds())
    else: log.info("Skipped.")

    if not sample_name:
        sample_name = os.path.basename(read1).split('_')[0].split('.')[0]

    log = rlog.getChild('report')
    start = datetime.datetime.now()

    log.critical(f"** {sample_name} is predicted to be {prediction}.**")

    note = []
    if checkpoint == 1:
        note.append("No read was mapped to the reference sequence database.")
    elif checkpoint == 2:
        note.append(f"{sample_name} is ipaH-.")
    elif checkpoint == 13:
        note.append("Shigella boydii serotype 13 is no longer considered a Shigella.")
    elif checkpoint == 31:
        note.append("No ipaH with sufficient coverage and accuracy was detected.")
    elif checkpoint == 32:
        note.append(f"{sample_name} is lacY+ or cadA+ but not one of the exception Shigella serotypes.")
    elif checkpoint == 41:
        note.append(f"No known wzx was detected. Either there was not enough coverage, or {sample_name} is a novel Shigella strain.")
    elif checkpoint == 42:
        note.append("Multiple wzx genes were detected. There's a potential contamination in the sample.")

    if ipaB >0:
        note.append("this strain is ipaB+, suggesting that it retains the virulent invasion plasmid.")

    note = ";".join(note)
    log.critical(note)
    # write hits table for individual sample
    
    if checkpoint == 1:
        pass
    else:
        output = f"{outdir}/{sample_name}-hits.tsv"
        log.critical(f"Writing hits to {output}")
        if "ipaH_c" in Maphits.Hit.tolist():
            List2.to_csv(output, sep="\t")
            List2.to_csv(sys.stdout, sep="\t")
        else:
            Maphits.to_csv(output, sep="\t")
            Maphits.to_csv(sys.stdout, sep="\t")

    serotype_output = f"{outdir}/{sample_name}.tsv"
    log.critical(f"Writing final serotype prediction to {serotype_output}")
    with open(serotype_output, "wt") as fh_out:
        ipab = ('-','+')[bool(ipaB)]
        results = f"sample\tprediction\tipaB\tnotes\n{sample_name}\t{prediction}\t{ipab}\t{note}"
        print(results)
        fh_out.write(f"{results}\n")

    lapse = datetime.datetime.now() - start
    log.critical(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

def main():
    ShigellaRef = os.path.abspath(j(dirname(__file__), 'resources', 'ShigellaRef5.fasta'))
    if not os.path.isfile(ShigellaRef):
        print(f"Error: reference sequence not found at {ShigellaRef}.", file=sys.stderr)
        exit(126)
    missing = []
    for dependency in ('minimap2', 'samtools', 'bcftools', ):
        try: 
            check_output([dependency, '--version'])
        except (CalledProcessError, FileNotFoundError) as e:
            missing.append((dependency, e))

    if missing:
        for dependency, error in missing:
            print(f"Dependency {dependency} errored or not found", file=sys.stderr)
        exit(126)

    import argparse

    parse = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
    parse.add_argument('--R1', metavar="FASTA", type=str,
                       help="Input FASTQ is R1 of paired-end reads")
    parse.add_argument('--R2', metavar="FASTA", type=str,
                       help="Input FASTQ is R2 of paired-end reads")
    parse.add_argument('--SE', metavar="FASTA", type=str,
                       help="Input FASTQ is contains single-end reads")
    parse.add_argument('--ont', action="store_true",
                        help="The input FASTQ file contains ONT reads")
    parse.add_argument('-n', '--name', dest='sample_name')
    parse.add_argument('-o', '--outdir', help='Where to write output files to', default="./")
    parse.add_argument('--verbose', '-v', action='count', dest='log_v', default=0)
    parse.add_argument('--version', action='version', version=f"ShigaTyper {__version__}")

    if len(sys.argv) == 1:
        parse.print_help()
        sys.exit(0)

    args = parse.parse_args()

    # Check inputs exist
    is_paired = False
    if args.R1 and args.R2 and args.SE:
        logging.error(f"Can only provide paired-end reads (--R1, --R2) or single-end reads (--SE), not both. Please fix and try again")
        sys.exit(1)
    elif args.R1 and args.R2:
        if not os.path.exists(args.R1):
            logging.error(f"Input R1 FASTQ ({args.R1}) does not exist, please verify and try again")
            sys.exit(1)
        if not os.path.exists(args.R2):
            logging.error(f"Input R2 FASTQ ({args.R2}) does not exist, please verify and try again")
            sys.exit(1)
        is_paired = True
    elif args.SE:
        if not os.path.exists(args.SE):
            logging.error(f"Input single-end FASTQ ({args.SE}) does not exist, please verify and try again")
            sys.exit(1)
    elif not args.R1 and not args.R2 and not args.SE:
        logging.error(f"A FASTQ file(s) must be provided, please correct and try again")
        sys.exit(1)

    for level, _ in zip((logging.ERROR, logging.INFO, logging.DEBUG), range(args.log_v + 1)):
        pass

    logging.basicConfig(
        level = level,
        stream=sys.stderr,
        format="[%(levelname)s::%(name)s] %(message)s"
    )

    try:
        tempdir = tempfile.mkdtemp()
        reads = []
        if is_paired:
            run([args.R1, args.R2], tempdir=tempdir, ont=False, sample_name=args.sample_name, outdir=args.outdir)
        else:
            run([args.SE], tempdir=tempdir, ont=args.ont, sample_name=args.sample_name, outdir=args.outdir)
    except CalledProcessError as e:
        rlog.error(e)
        if e.stderr:
            rlog.error(e.stderr)
        exit(e.returncode)
    finally:
        pass
        shutil.rmtree(tempdir)

if __name__ == "__main__":
    main()
