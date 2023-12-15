# Cell-Flex-multi (singleplex)
Here we are going to analyze a new 10x genomics in the lab called cellflex, Fixation with paraformaldehyde (PFA) allows samples to be collected, shipped to a central location, and analyzed without sacrificing integrity or data quality, creating new possibilities for sample accessibility, throughput, and batched analysis. you can learn more about this new kit from here: https://www.10xgenomics.com/products/single-cell-gene-expression-flex

We are going to analyze SCGTEST_51 which is the first test of this kit, it is 2 samples and we are going to follow this toutorial from 10X: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi-frp#probe-set-files

To work with this we need to get lims information, which is all thing about the project we are analyze 

# 1- Get LIMS info

### The script 
```
#!/bin/bash

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project
/home/groups/singlecell/mabdalfttah/projects/scripts/limsq.py -sp $1 | sed 's/;/                                                                                                                                                             \t/g' > "lims_info_"$1".txt"

echo "Created LIMS information file: lims_info.txt"
```
### How to run

**Note**: make sure you write the righ path of a script called limsq.py which is the main script to get the lins info

```
./1-lims.sh SCGTEST_49
```
Now we have lims_info_SCGTEST_51.txt with all of information of the samples and project

Since some fastq files/samples doesn't pass the filters of QC we need to remove them from the lims_info_SCGTEST_51.txt file 

```
awk -F'\t' -v column="LanePassFail" 'BEGIN {OFS=FS} NR==1 {for (i=1; i<=NF; i++) if ($i == column) col=i} $col != "fail"' lims_info_SCGTEST_51.txt > tmp_file && mv tmp_file lims_info_SCGTEST_51.txt
```

We have filterd lims file, we are ready for the next step 

# 2- Get FASTQs Path

We need to know where is the FASTQs path of our project, to have this we need to run this script, which is use the lims info file and write the path of the FASTQs, usually the FASTQs are in this path **/scratch/project/production/fastq** if they are not you need to change this path 

### The script 
```
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
```
### How to run 

**NOTE**: First we need to activate any conda env with python:
```
source ~/.bashrc
conda activate sc_py
```
Run the script 
```
python 2-write_fastq_paths.py --subproject SCGTEST_51 --info_file lims_info_SCGTEST_51.txt
```

# 3- Create a Metadata File
This step is about create a metadata file from lims infor file, it is just information about our samples and which library they they are in, it is a different format of LIMS with just the important information

**This step should be in R**

```
subproject = "SCGTEST_52"
data <- read.table(glue::glue("Y:/projects/{subproject}/lims_info_{subproject}.txt"), header = T, sep = "\t")

metadata = data.frame(subproject = unique(data$subproject),
                      gem_id = unique(data$SampleName),
                      library_id = unique(data$id),
                      library = unique(data$library),
                      type = "not_hashed",
                      donor_id = unique(data$SampleName),
                      flowcell = unique(data$flowcell),
                      lane = unique(data$lane),
                      index = unique(data$index))

write.csv(x = metadata, glue::glue("Y:/projects/{subproject}/metadata.csv"), row.names = F)
```
Now we have SCGTEST_51.csv metadata file, let's go for the next step

# 4- Create a jobs directories and copy FASTQs to them

In this script we create a directory for each samples (We know from metadata) and copy the FASTQs files to this directory, wiht this we will have the FASTQs of each sample in seperate diretory to pass this directory when we run cellranger (it run per sample seperatly)

### The script 
```
# This script initializes the filesystem of this project:
# It creates a "jobs" folder which contains as many subdirectories as samples it has
# For each sample directory, it creates the following files/folders:
# 1. fastq: dir with the symlinks pointing to the fastq files
# 2. log: dir which contains standard error and output of cellranger
# 3. (sample_id).cmd: job script to compute the features-barcode matrix using cellranger


# Import required packages
import numpy as np
import pandas as pd
import os
import argparse
import subprocess
import re
import sys
import config_vars as cfg
from utils import *


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--gem_id",
                    dest = "gem_id",
                    action = "store",
                    default = None,
                    help = "Gel Beads in Emulsion id")
parser.add_argument("--verbose",
                    dest = "verbose",
                    action = "store_true",
                    default = False,
                    help = "Print log in standard error")
parser.add_argument("--metadata",
                    dest = "metadata",
                    action = "store",
                    default = None,
                    help = "Metadata csv file for the tonsil atlas project")
parser.add_argument("--fastq_paths",
                    dest = "fastq_paths",
                    action = "store",
                    default = None,
                    help = "File that contains the paths of the fastqs for the subproject libraries")


def create_fastq_symlink_nh(gem_id, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation

    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      fastq_path_df: pandas dataframe with the fastq paths for that gem_id
      symlink_path: string specifying where to create the symlinks

    Returns:
      None
    """
    pair_ids = np.unique(fastq_path_df["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            fastq_path = pair_df.loc[j, "fastq_path"]
            lane = str(i + 1)
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path, gem_id, lane, read)])

options = parser.parse_args()
subproject = options.subproject
gem_id = options.gem_id
metadata_path = options.metadata
fastq_paths = options.fastq_paths


# Read data
project_dir = "/home/groups/singlecell/mabdalfttah/projects/{}".format(subproject)
fastq_path_df = pd.read_csv(fastq_paths, sep = "\t", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
if options.verbose:
    sys.stderr.write("Files read successfully!\n")


# For each sample, create directories and jobscript
if not os.path.exists("{}/jobs".format(project_dir)):
    os.mkdir("{}/jobs".format(project_dir))
filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]


# Create directories
subproject_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(subproject_dir)
log_dir = "{}/log".format(subproject_dir)
for direct in [subproject_dir, fastq_dir, log_dir]:
    if not os.path.exists(direct):
        os.mkdir(direct)


# Define variables and subset dataframes
library_id = metadata_df.loc[filt, "library_id"]
fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_id), :]
type = metadata_df["type"]
type = type.values[0]

# Create symmlinks to fastq files
create_fastq_symlink_nh(gem_id, fastq_sub_df, fastq_dir)
```

### How to run
```
python 3-copy_fastqs  --subproject SCGTEST_51 --fastq_paths fastq_paths.tab --metadata SCGTEST_51.csv --gem_id CNAG_61_CellFlex_A
python 3-copy_fastqs  --subproject SCGTEST_51 --fastq_paths fastq_paths.tab --metadata SCGTEST_51.csv --gem_id CNAG_61_CellFlex_B
```

# Create Config csv file 
All pevoius steps was internal steps to make the FASTQs ready to analyze, but from here the steps are realted to cellranger irself, the first thing we need to do for cellranger is to create a config file.
The multi config CSV contains both the library definitions and experimental design variables. The required sections differ slightly for analysis with single-sample ("singleplex") vs. multiplexed configurations. It is composed of up to four sections for FRP data:

-  The [gene-expression] section has two columns that specify parameters relevant to analysis of gene expression data, such as reference genome and cell-calling parameters that apply to the whole library in the case of singleplex Fixed RNA Profiling.
- The [feature] section has two columns that specify parameters relevant to analysis of Feature Barcode libraries; required when analyzing singleplex Fixed RNA Profiling data with Cell Surface Protein ("Antibody Capture" library).
- The [libraries] section has three required columns that specify where the input FASTQ files may be found.
- The [samples] section has two required columns that specify sample information for Fixed RNA Gene Expression with multiple samples or a single sample with multiple Probe Barcodes; it is not valid for **singleplex** Fixed RNA Profiling analysis.

**NOTE**: Our analysis here is for singleplex GEX 

You can know more about this from here https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi-frp#probe-set-files

### The script 

```{}
import csv
import sys
import os

# Check if four command-line arguments are provided
if len(sys.argv) != 5:
    print("Usage: python create_config_csv.py <cells> <Sample> <FASTQs> <subproject>")
    sys.exit(1)

# Accept user inputs
cells = sys.argv[1]
sample = sys.argv[2]
fastqs = sys.argv[3]
subproject = sys.argv[4]

# Define the data with user inputs
data = [
    ["[gene-expression]",],
    ["reference","/scratch/groups/singlecell/data/reference/refdata-gex-GRCh38-2020-A",],
    ["probe-set","/scratch/groups/singlecell/software/cellranger/7.1.0/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv",],
    ["no-bam", "true",],
    ["expect-cells", cells,],
    ["[libraries]",],
    ["fastq_id", "fastqs", "feature_types"],
    [sample, fastqs, "Gene Expression"],
]

# Define the output file path
output_file_path = f"/home/groups/singlecell/mabdalfttah/projects/{subproject}/jobs/{sample}/config.csv"

# Ensure the directory structure exists
os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

# Write the data to the specified CSV file path
with open(output_file_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerows(data)

print(f"CSV file '{output_file_path}' created successfully.")
```

### How to run 
There is four variables in this script:

1- The number of expected cells (Ask wet lab)

2- path of FASTQs (you should know this when you run copy FASTQs script)

3- subproject and  sample (gem_id) [it is use to save the config.csv file in this path /home/groups/singlecell/mabdalfttah/projects/{subproject}/jobs/{sample}/config.csv]


```{}
python 4-create_config_csv.py 8000 CNAG_61_CellFlex_A /home/groups/singlecell/mabdalfttah/projects/SCGTEST_51/jobs/CNAG_61_CellFlex_A/fastq  SCGTEST_51
python 4-create_config_csv.py 8000 CNAG_61_CellFlex_B /home/groups/singlecell/mabdalfttah/projects/SCGTEST_51/jobs/CNAG_61_CellFlex_B/fastq  SCGTEST_51
```

# Create a job file
Finaly last but least we will create a job script to run cellranger on the cluster, this job script contain an inportant information for the cluster and the cellranger command.

**Note**: since this kit in new kit, so the new cellranger versions of cellranger used to analyze this data, here we use cellranger 7.1.0, and it is installed in **/scratch/groups/singlecell/software/cellranger/7.1.0/cellranger**

**NOTE**: we use cellranger multi not count even for singleplex fixed Gene-expression 

```{}
#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <sample> <subproject>"
    exit 1
fi

# Accept user inputs
sample="$1"
subproject="$2"

# Define the output file path
output_file_path="/home/groups/singlecell/mabdalfttah/projects/${subproject}/jobs/${sample}/${sample}.cmd"

# Create the content for the .cmd file
cat > "$output_file_path" <<EOF
#!/bin/bash
#SBATCH --job-name=$sample

#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=mohamed.abdalfttah@cnag.crg.eu

#SBATCH --output=%x.slurm.%J.out        # define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err

#SBATCH --time=11:00:00 # set a maximum time that the job will take HH:MM:SS (process will be terminated after this is reached)

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --partition=genD
#SBATCH --mem=32G

echo [\`date "+%Y-%m-%d %T"\`] started job on \$HOSTNAME

export TENX_IGNORE_DEPRECATED_OS=1
export HDF5_USE_FILE_LOCKING=FALSE

/scratch/groups/singlecell/software/cellranger/7.1.0/cellranger multi --id $sample --csv config.csv  --jobmode /scratch/groups/singlecell/software/cellranger/6.1.1/external/martian/jobmanagers/new_cluster/slurm.template;
EOF

echo "Created .cmd file: $output_file_path"
```

### How to run 

```{}
chmod +x create_cmd_file.sh
./5-create_job_script.sh CNAG_61_CellFlex_A SCGTEST_51
./5-create_job_script.sh CNAG_61_CellFlex_B SCGTEST_51
```

Finally we submit the jobs and going home 

```{}
cd jobs/CNAG_61_CellFlex_A
sbatch CNAG_61_CellFlex_A.cmd
cd jobs/CNAG_61_CellFlex_B
sbatch CNAG_61_CellFlex_B.cmd
```

When cellranger finish you will found a directory with the name of the sample contain all of the result from the run 
