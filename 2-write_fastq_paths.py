#!/usr/bin/env python

# Writes fastq path by arranging proper flowcell, lane, index and read for a set of libraries

# Load packages
import numpy as np
import pandas as pd
import os
import argparse


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to transfer feature-barcodes matrices from cluster to lcal")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--info_file",
                    dest = "info_file",
                    action = "store",
                    default = None,
                    help = "Tab-delimited file with the information of Illumina sequence of libraries for that subproject")


options = parser.parse_args()
subproject = options.subproject
info_file = options.info_file

# Read file
lims = pd.read_csv(info_file, sep = "\t", header = 0)

# Assemble fastq paths combining flowcell, lane and index
fastq_path = "/scratch/project/production/fastq"
fastq_path_list_r1 = []
fastq_path_list_r2 = []
for idx in lims.index:
    fc = lims.loc[idx, "flowcell"]
    lane = lims.loc[idx, "lane"]
    index = lims.loc[idx, "index"]
    fastq_path_r1 = "{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_r2 = "{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_list_r1.append(fastq_path_r1)
    fastq_path_list_r2.append(fastq_path_r2)
library_id_l = list(lims["id"].append(lims["id"]))
p_l = "P" * len(fastq_path_list_r1)
indx_l = list(range(1, len(fastq_path_list_r1) + 1))
pair_id = [p_l[x] + str(indx_l[x]) for x in range(len(indx_l))]
fastq_path_list_r1.extend(fastq_path_list_r2)
pair_id.extend(pair_id)
fastq_path_l = fastq_path_list_r1
read_l = (["R1"] * lims.shape[0]) + (["R2"] * lims.shape[0])
fastq_dict = {"library_id":library_id_l, "fastq_path":fastq_path_l, "read":read_l, "pair_id":pair_id}
fastq_df = pd.DataFrame(fastq_dict)

fastq_df.to_csv("fastq_paths.tab".format(subproject), header = True, index = False, sep="\t")
