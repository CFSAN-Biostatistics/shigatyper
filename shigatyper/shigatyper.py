#!/usr/bin/env python3
# coding: utf-8

version = "1.0.5"

usage = f"""ShigaTyper v. {version}, 2019

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
import glob
import logging
import itertools
#import numpy as np
import pandas as pd
import shutil
import sys
import tempfile

from subprocess import check_output, CalledProcessError
from os.path import join as j, dirname
from functools import partial




#from IPython.display import display, HTML

# start = datetime.datetime.now()
# print("read1: ", read1)
# print("read2: ", read2)
# print("name of sample: ", Sample)
# #create a directory for all files generated from analyzing this sample
# if not os.path.isdir(Sample):
#     os.makedirs(Sample)
# timetrack = []
# lapse = datetime.datetime.now() - start; timetrack.append(lapse.total_seconds())


# # ## 2. QC Characterization of the two reads files by fastp

# # ### 2.1. summary of quality attributes of the two fastq files (read1 and read2), no filtering:

# # In[ ]:


# start = datetime.datetime.now()
# today = datetime.datetime.now().strftime("%Y%m%d")
# if fastp == 0:
#     Report_html = Sample + "/" + Sample + "_fastp_" + today + ".html"
#     Report_json = Sample + "/" + Sample + "_fastp_" + today + ".json"
#     fastplog = Sample + "/" + Sample + "fastp.log"
#     get_ipython().system('fastp -i $read1 -I $read2 -Q -h $Report_html -j $Report_json >> $fastplog')
# else: print("Quality insepction skipped.")
# lapse = datetime.datetime.now() - start; timetrack.append(lapse.total_seconds())


# # In[ ]:


# start = datetime.datetime.now()
# if fastp == 0: 
#     import json
#     data = json.load(open(Report_json))

#     attr = ["Number of reads", "Number of bases", "Q20 bases", "Q30 bases", "Average read length"]
#     read1sum = [data["read1_before_filtering"].get("total_reads"), data["read1_before_filtering"].get("total_bases"),
#                 data["read1_before_filtering"].get("q20_bases"), data["read1_before_filtering"].get("q30_bases")]
#     read1q20 = str(read1sum[2]) + " (" + str(round(100*read1sum[2]/read1sum[1],2)) +"%)"
#     read1q30 = str(read1sum[3]) + " (" + str(round(100*read1sum[3]/read1sum[1],2)) +"%)"
#     col1 = read1sum[0:2] + [read1q20, read1q30] + [round(read1sum[1]/read1sum[0],2)]
#     read2sum = [data["read2_before_filtering"].get("total_reads"), data["read2_before_filtering"].get("total_bases"),
#                 data["read2_before_filtering"].get("q20_bases"), data["read2_before_filtering"].get("q30_bases")]
#     read2q20 = str(read2sum[2]) + " (" + str(round(100*read2sum[2]/read2sum[1],2)) +"%)"
#     read2q30 = str(read2sum[3]) + " (" + str(round(100*read2sum[3]/read2sum[1],2)) +"%)"
#     col2 = [read2sum[0], read2sum[1], read2q20, read2q30] + [round(read2sum[1]/read2sum[0],2)]
#     totalsum = [data["summary"]["before_filtering"].get("total_reads"), data["summary"]["before_filtering"].get("total_bases"),
#                 data["summary"]["before_filtering"].get("q20_bases"), data["summary"]["before_filtering"].get("q30_bases")]
#     totalq20 = str(totalsum[2]) + " (" + str(round(100*totalsum[2]/totalsum[1],2)) +"%)"
#     totalq30 = str(totalsum[3]) + " (" + str(round(100*totalsum[3]/totalsum[1],2)) +"%)"
#     col3 = totalsum[0:2] + [totalq20, totalq30] + [round(totalsum[1]/totalsum[0],2)]

#     df = pd.DataFrame({' ': attr,'read1': col1,'read2':col2, 'Total': col3})
#     df = df[[' ', 'read1', 'read2', 'Total']]
    
#     display(HTML(df.to_html(index=False)))
# lapse = datetime.datetime.now() - start; timetrack.append(lapse.total_seconds())


# # ### 2.2. Visualization of base quality by type and position

# # In[ ]:


# start = datetime.datetime.now()
# get_ipython().run_line_magic('matplotlib', 'inline')
# if fastp == 0:
#     read1_after = data["read1_before_filtering"]["quality_curves"]
#     read1A = read1_after["A"]; read1G = read1_after["G"]
#     read1T = read1_after["T"]; read1C = read1_after["C"]

#     read2_after = data["read2_before_filtering"]["quality_curves"]
#     read2A = read2_after["A"]; read2G = read2_after["G"]
#     read2T = read2_after["T"]; read2C = read2_after["C"]

#     import matplotlib.pyplot as plt
#     fig, (plt1, plt2) = plt.subplots(nrows =1, ncols =2, figsize = (6*2, 5*1))

#     plt1.plot(range(len(read1A)), read1A, color = "green")
#     plt1.plot(range(len(read1G)), read1G, color = "black")
#     plt1.plot(range(len(read1T)), read1T, color = "red")
#     plt1.plot(range(len(read1C)), read1C, color = "blue")
#     plt1.set_ylim((15, 43)); plt1.set_title("Read1: Base Quality")
#     plt1.legend(("read1_A", "read1_G", "read1_T", "read1_C"), loc = "lower left")
#     plt1.set_xlabel("Position"); plt1.set_ylabel("Average Quality Score")

#     plt2.plot(range(len(read2A)), read2A, color = "green")
#     plt2.plot(range(len(read2G)), read2G, color = "black")
#     plt2.plot(range(len(read2T)), read2T, color = "red")
#     plt2.plot(range(len(read2C)), read2C, color = "blue")
#     plt2.set_ylim((15, 43)); plt2.set_title("Read2: Base Quality")
#     plt2.legend(("read2_A", "read2_G", "read2_T", "read2_C"), loc = "lower left")
#     plt2.set_xlabel("Position")
# else: print("Not plotted.")
# lapse = datetime.datetime.now() - start; timetrack.append(lapse.total_seconds())


# # ### 2.3. Average depth of coverage

# # In[ ]:


# start = datetime.datetime.now()
# if fastp == 0:
#     print("Depth of coverage (Assuming a genome size of ~5 Mbp): ", round(totalsum[1]/5e6, 1), " fold")
# else: print("Not calculated.")
# lapse = datetime.datetime.now() - start; timetrack.append(lapse.total_seconds())

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

def run(read1, read2, tempdir, sample_name = '', threshold=20, rlog=rlog, *args, **kwargs):
    # ## 3. Map Filtered reads to Reference sequence database (ShigellaRef5)

    rlog.critical(f" .v {version}")

    rlog = rlog.getChild('run')
    
    timetrack = []
    #sub = partial(check_output, shell=True)
    def sub(cmd, log=rlog):
        log = log.getChild(cmd.split()[0]) #log with name of subprocess command
        result = str(check_output(cmd, shell=True, encoding='UTF-8')).strip()
        log.info(cmd)
        for row in result.splitlines():
            log.debug("{}".format(row.replace(' ', '\t')))
        return result

    read1 = os.path.abspath(read1)
    read2 = os.path.abspath(read2)

    # In[ ]:

    
    log = rlog.getChild('index')
    start = datetime.datetime.now()
    # give an error message if the reference sequence database is not there
    ShigellaRef = os.path.abspath(j(dirname(__file__), 'resources', 'ShigellaRef5.fasta'))
    if not os.path.isfile(ShigellaRef):
        log.error(f"Error: reference sequence not found at {ShigellaRef}")
        exit()
    # indexing reference sequence database when necessary
    #dir_path = os.path.dirname(os.path.realpath(ShigellaRef)) # absolute path of reference sequence directory
    #rel_dir = os.path.relpath(dir_path, os.getcwd()) # relative path of reference sequence directory
    mmi_index = j(tempdir, "ShigellaRef5.mmi")
    if not os.path.isfile(mmi_index):
        log.warning(f"building Reference sequence index in {mmi_index}")
        sub(f'minimap2 -d {mmi_index} {ShigellaRef}', log)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    # In[ ]:

    log = rlog.getChild('map')
    start = datetime.datetime.now()
    # map the fastq.gz files to reference sequence database
    outputbam = j(tempdir, "_Shigella5.bam")
    #maplog = j(tempdir, "_Shigella5_minimap2.log")
    sub(f'minimap2 -ax sr {mmi_index} {read1} {read2} | samtools view -F 0x04 -b | samtools sort -o {outputbam}', log)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())


    # ### Checkpoint 1 
    # If there is no read mapped to the reference sequence database, discontinue the analysis.

    # In[ ]:

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


    # ## 3.1. Examine sequence hits identified by bowtie2

    # In[ ]:

    log = rlog.getChild('hits')
    start = datetime.datetime.now()
    #check what sequences were hits and how many reads were mapped to each of the hits
    if checkpoint == 0:
        Hits = sub(f'samtools view {outputbam} | cut -f3 | uniq -c', log)
        hits = []
        Nreads = []
        # for hit in Hits:
        #     hits.append(hit[8:]); Nreads.append(int(hit[:7]))
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

    # In[ ]:

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

    # In[ ]:

    log = rlog.getChild('accuracy')
    start = datetime.datetime.now()
    if checkpoint == 0:
        #print("........................")
        # find reference sequence length for calculation of % coverage
        Gene_length = list()
        RefDic = dict()
        Reflines = sub(f'samtools view -H {outputbam} | grep "SN"', log)
        # for Refline in Reflines:
        #     line = Refline[(Refline.find("SN:")+3):]
        #     words = line.split("LN:")
        #     RefDic[words[0].strip()] = int(words[1])
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


    # In[ ]:

    log = rlog.getChild('pileup')
    start = datetime.datetime.now()
    if checkpoint == 0:
        #log.info("Analysis in progress........ \n")
        # index the bam file for mpileup
        sub(f'samtools index {outputbam}', log)
        # index reference sequence for mpileup if the index is not already there
        #if os.path.isfile("../../references/ShigellaRef5.fasta.fai") == False:
        #    !samtools faidx $ShigellaRef

        #outputmpileup = Sample + "/" + Sample + "_Shigella5.mpileup"
        outputmpileup = j(tempdir, "Shigella5.mpileup")
        sub(f'samtools mpileup -C50 -q 20 -Q 20 -f {ShigellaRef} {outputbam} -o {outputmpileup}', log)
        # can I pipe it so that I won't have a mpileup file left in the disk?
        lapse = datetime.datetime.now() - start
        log.info(f"Complete in {readable(lapse)}.")
        timetrack.append(lapse.total_seconds())
    else: log.info("skipped.")



    # In[ ]:

    log = rlog.getChild('hits')
    start = datetime.datetime.now()
    if checkpoint == 0:
        covSummary = sub(f"cat {outputmpileup} | awk '$4 > 0' | cut -f1 | uniq -c", log)
        hits = []; bpCovered = []
        # for hit in covSummary:
        #     hits.append(hit[8:]); bpCovered.append(int(hit[:7]))
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


    # In[ ]:

    log = rlog.getChild('variants')
    start = datetime.datetime.now()
    if checkpoint == 0:
        variants = sub(f'samtools mpileup -C50 -q 20 -Q 20 -f {ShigellaRef} -g {outputbam} | bcftools     call -m | cat | grep -v "^#" | grep PL | cut -f1 | uniq -c', log)
        #VARS = VARS[3:]
        hits = []; Nvar = []
        # for hit in VARS:
        #     hits.append(hit[8:]); Nvar.append(int(hit[:7]))
        #     Variants = pd.DataFrame({'Hit': hits, 'Number of variants': Nvar})
        # 7 ipaH_c
        # 17 ipaB
        for nvar, hit in pairwise(variants.split()):
            hits.append(hit)
            Nvar.append(int(nvar))
        Variants = pd.DataFrame({'Hit':hits, 'Number of variants':Nvar})

    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())


    # In[ ]:

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


    # In[ ]:

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

    # In[ ]:

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

    # In[ ]:

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
                # np.delete(Hits, wzxindex(Hits)).tolist() #import numpy as np # I forgot why I have this line
                checkpoint = 42
                prediction = "No prediction (multiple wzx)"
    else: log.info("Skipped.")
                
    if checkpoint == 0: log.info("Checkpoint 4 ..... passed")
    elif checkpoint == 4: log.error("Checkpoint 4 failed!")
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())


    # ## 4. Shigella serotype prediction

    # In[ ]:

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
                    prediction = "Shigalla dysenteriae serotype 1, rfp- (phenotypically negative)"
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
                    a = [int(i) for i in a]; average_a = sum(a)/len(a)
                    if a == 0:
                        prediction = "Shigella boydii serotype 6 or 10"
                    else:
                        b = sub(f"cat {outputmpileup} | awk '$1==\"wbaM\" && $2 >251 && $2 <254' | cut -f4")
                        b = [int(i) for i in b]; average_b = sum(b)/len(b)
                        junction_ratio = average_b/average_a
                        c = sub(f'samtools depth -r wbaM {outputbam} | cut -f3')
                        c = [int(i) for i in c]; average_c = sum(c)/len(c)
                        d = sub(f'cat {outputmpileup} | awk \'$1=="wbaM"\' | cut -f4')
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
                        Hits.remove("Sf_wzx")
                        try: Hits.remove("Sf_wzy")
                        except: Hits = Hits
                        if len(Hits) == 0:
                            prediction = "Shigella flexneri serotype Y"
                        else:
                            SfDic = {"Shigella flexneri Yv": ["Xv"], "Shigella flexneri serotype 1a": ["gtrI"], 
                                "Shigella flexneri serotype 1b": ["gtrI", "Oac1b"], "Shigella flexneri serotype 2a": 
                                ["gtrII"], "Shigella flexneri serotype 2b": ["gtrII", "gtrX"], "Shigella flexneri serotype 3a":
                                ["gtrX","Oac"], "Shigella flexneri serotype 3b": ["Oac"], "Shigella flexneri serotype 4a": 
                                ["gtrIV"], "Shigella flexneri serotype 4av": ["gtrIV", "Xv"], "Shigella flexneri serotype 4b":
                                ["gtrIV", "Oac"], "Shigella flexneri serotype 5a": (["gtrV", "Oac"], ['gtrV']),
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



    # In[ ]:
    if not sample_name:
        sample_name = os.path.basename(read1).split('_')[0].split('.')[0]

    log = rlog.getChild('report')
    start = datetime.datetime.now()
    # import papermill as pm
    # pm.record("prediction", prediction)

    #from IPython.display import Markdown, display
    log.critical(f"** {sample_name} is predicted to be {prediction}.**")
    #display(Markdown(finaloutput))
    if checkpoint == 1:
        log.critical("No read was mapped to the reference sequence database.")
    elif checkpoint == 2:
        log.critical("{sample_name} is ipaH-.")
    elif checkpoint == 13:
        log.critical("Shigella boydii serotype 13 is no longer considered a Shigella.")
    elif checkpoint == 31:
        log.critical("No ipaH with sufficient coverage and accuracy was detected.")
    elif checkpoint == 32:
        log.critical("{sample_name} is lacY+ or cadA+ but not one of the exception Shigella serotypes.")
    elif checkpoint == 41:
        log.critical(f"No known wzx was detected. Either there was not enough coverage, or {sample_name} is a novel Shigella strain.")
    elif checkpoint == 42:
        log.critical("Multiple wzx genes were detected. There's a potential contamination in the sample.")


    if ipaB >0:
        log.critical("this strain is ipaB+, suggesting that it retains the virulent invasion plasmid.")

    # if checkpoint == 1:
    #     pass
    # else:
    #     log.critical("Please consult the table below for further information:")
    #     from IPython.display import display, HTML
    #     if "ipaH_c" in Maphits.Hit.tolist():
    #         rowindex = []
    #         for i in range(List2.shape[0]):
    #             if List2.iloc[i]['% covered'] > threshold:
    #                 rowindex.append(i)
    #         def color(x):
    #             df = x.copy()
    #             df.loc[:,:] = ""
    #             df.loc[rowindex, ] = 'color: blue'
    #             return df
    #         List2_blue = List2.style.apply(color, axis = None)
    #         display(HTML(List2_blue.render(index=False)))
    #         log.critical(f"Note: colored in blue are gene hits that passed threshold length coverage. ({threshold}% )")
    #     else: display(HTML(Maphits.to_html(index=False)))
            
    #now = datetime.datetime.now().strftime("%Y-%m-%d %H-%M")
    #log.info("\nDate and time of analysis: ", now)
    lapse = datetime.datetime.now() - start
    log.info(f"Complete in {readable(lapse)}.")
    timetrack.append(lapse.total_seconds())

    wtr = csv.writer(sys.stdout, delimiter='\t', dialect='excel')
    wtr.writerow(('sample', 'prediction', 'ipaB'))
    wtr.writerow((sample_name, prediction, ('-','+')[bool(ipaB)]))


    # In[ ]:


    # # removing files generated from fastp, minimap2, and samtools:
    # start = datetime.datetime.now()
    # filedir = Sample + "/*.*"
    # files = glob.glob(filedir)
    # for file in files: os.remove(file)
    # os.rmdir(Sample)
    # lapse = datetime.datetime.now() - start
    # timetrack.append(lapse.total_seconds())


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
    parse.add_argument('read1')
    parse.add_argument('read2')
    parse.add_argument('-n', '--name', dest='sample_name')
    parse.add_argument('--verbose', '-v', action='count', dest='log_v', default=0)
    parse.add_argument('--version', action='version', version=f"ShigaTyper v. {version}")

    


    args = parse.parse_args()


    for level, _ in zip((logging.ERROR, logging.INFO, logging.DEBUG), range(args.log_v + 1)):
        pass

    logging.basicConfig(
        level = level,
        stream=sys.stderr,
        #format="%(asctime)s %(levelname)s:%(name)s:\t%(message)s"
        format="[%(levelname)s::%(name)s] %(message)s"
    )


    try:
        tempdir = tempfile.mkdtemp()
        run(tempdir=tempdir, **vars(args))
    except CalledProcessError as e:
        #print(e.cmd, file=sys.stderr)
        rlog.error(e)
        if e.stderr:
            rlog.error(e.stderr)
        exit(e.returncode)
    finally:
        pass
        shutil.rmtree(tempdir)

if __name__ == "__main__":
    main()

# save output for time
#if fastp== 0:
#    import csv
#    timeoutput = Sample + "_" + now[:now.find(" ")] + "_time.csv"
#    with open(timeoutput, 'w') as f:
#        writer = csv.writer(f)
#        writer.writerow([Sample] + timetrack)



