# Python functions used in different scripts

import numpy as np
import pandas as pd
import subprocess
import os
import config_vars as cfg

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



def make_cellranger_nh(gem_id, jobscript_path, fastq_path, expected_cells):
    """Creates a cellranger script for a non-hashed GEM well 
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      fastq_path: path to the fastq files
      expected_cells: expected number of high-quality cells in this experiment

    Returns:
      None 
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash
#SBATCH --job-name={}    

#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=mohamed.abdalfttah@cnag.crg.eu

#SBATCH --output=%x.slurm.%J.out 	# define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err

#SBATCH --time=23:00:00	# set a maximum time that the job will take HH:MM:SS (process will be terminated after this is reached)


#SBATCH --cpus-per-task=1
#SBATCH --partition=genB,main
#SBATCH --ntasks=1
#SBATCH --qos=normal

export TENX_IGNORE_DEPRECATED_OS=1 
export HDF5_USE_FILE_LOCKING=FALSE 

    
{} count --fastqs {} --id {} --transcriptome {} --jobmode /scratch/groups/hheyn/software/cellranger/6.1.1/external/martian/jobmanagers/slurm.template;
""".format(gem_id, cfg.cellranger_path, fastq_path, gem_id, cfg.reference_path)
    job_script_file.write(job_script)
    job_script_file.close()


def create_fastq_symlink_h(gem_id, library, lib_type, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation for cell-hashed samples.
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      library: library id
      lib_type: type of library (cDNA or HTO)
      fastq_path_df: pandas dataframe with the fastq paths for that gem_id
      symlink_path: string specifying where to create the symlinks (fastq/HTO or fastq/cDNA)
    
    Returns:
      None 
    """
    fastq_path_sub = fastq_path_df.loc[fastq_path_df["library_id"] == library, :]
    pair_ids = np.unique(fastq_path_sub["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            lane = str(i + 1)
            symlink_path_lane = "{}/lane{}".format(symlink_path, lane)
            if not os.path.exists(symlink_path_lane):
                os.mkdir(symlink_path_lane)
            fastq_path = pair_df.loc[j, "fastq_path"]
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            if lib_type == "hashed_hto":
                gem_id_sp = "{}_HTO".format(gem_id)
            elif lib_type == "hashed_cdna":
                gem_id_sp = gem_id
            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path_lane, gem_id_sp, lane, read)])


def write_libraries_csv(gem_id, gem_id_path):
    """Creates the file "libraries.csv" which is required by cellranger in feature-barcoding analysis.
    
    Args:
     gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
     fastq_path: path to GEM-specific directory
     gem_id_path: absolute path to gem_id-specific directory.
    
    Returns:
      None 
    """
    lib_csv = open("{}/libraries.csv".format(gem_id_path), "w")
    lib_csv.write("fastqs,sample,library_type")
    fastq_dirs = os.listdir("{}/fastq".format(gem_id_path))
    for d in fastq_dirs:
        if d == "HTO":
            gem_id_sp = "{}_HTO".format(gem_id)
            lib_type = "Antibody Capture"
        elif d == "cDNA":
            gem_id_sp = gem_id
            lib_type = "Gene Expression"
        fastq_sub_dirs = os.listdir("{}/fastq/{}".format(gem_id_path, d))
        for sub_d in fastq_sub_dirs:
            sub_d_abs_path = "{}/fastq/{}/{}".format(gem_id_path, d, sub_d)
            output_line = "\n{},{},{}".format(sub_d_abs_path, gem_id_sp, lib_type)    
            lib_csv.write(output_line)
    lib_csv.close()


def make_cellranger_h(gem_id, jobscript_path, expected_cells):
    """Creates a cellranger script for a hashed GEM well 
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      expected_cells: expected number of high-quality cells in this experiment

    Returns:
      None 
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash

#SBATCH --job-name={}
#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=will.blevins@cnag.crg.eu
#SBATCH --output=%x.slurm.%J.out        # define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err
#SBATCH --time=23:00:00 # set a maximum time that the job will take HH:MM:SS (process will be terminated after this is >
#SBATCH --cpus-per-task=1
#SBATCH --partition=genB,main
#SBATCH --ntasks=1
#SBATCH --qos=normal

export TENX_IGNORE_DEPRECATED_OS=1 
export HDF5_USE_FILE_LOCKING=FALSE 

{} count --libraries libraries.csv --feature-ref feature_reference.csv --id {} --transcriptome {} --jobmode /scratch/groups/hheyn/software/cellranger/6.1.1/external/martian/jobmanagers/slurm.template;
""".format(gem_id, cfg.cellranger_path, gem_id, cfg.reference_path)
    job_script_file.write(job_script)
    job_script_file.close()


def make_cellranger_cite_seq(gem_id, jobscript_path, expected_cells):
    """Creates a cellranger script for a CITE-seq experiment 
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      expected_cells: expected number of high-quality cells in this experiment

    Returns:
      None 
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash

#SBATCH --job-name={}    
#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=will.blevins@cnag.crg.eu
#SBATCH --output=%x.slurm.%J.out        # define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err
#SBATCH --time=23:00:00 # set a maximum time that the job will take HH:MM:SS (process will be terminated after this is >
#SBATCH --cpus-per-task=1
#SBATCH --partition=genB,main
#SBATCH --ntasks=1
#SBATCH --qos=normal

export TENX_IGNORE_DEPRECATED_OS=1 
export HDF5_USE_FILE_LOCKING=FALSE 


{} count --libraries libraries.csv --feature-ref feature_reference.csv --id {} --transcriptome {} --jobmode /scratch/groups/hheyn/software/cellranger/6.1.1/external/martian/jobmanagers/slurm.template;
""".format(gem_id, cfg.cellranger_path, gem_id, cfg.reference_path)
    job_script_file.write(job_script)
    job_script_file.close()


def make_cellranger_vdj(gem_id, jobscript_path, fastq_path):
    """Creates a cellranger script for a single-cell VDJ experiment 
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      expected_cells: expected number of high-quality cells in this experiment

    Returns:
      None 
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash

#SBATCH --job-name={}    
#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=will.blevins@cnag.crg.eu
#SBATCH --output=%x.slurm.%J.out        # define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err
#SBATCH --time=23:00:00 # set a maximum time that the job will take HH:MM:SS (process will be terminated after this is >
#SBATCH --cpus-per-task=1
#SBATCH --partition=genB,main
#SBATCH --ntasks=1
#SBATCH --qos=normal

export TENX_IGNORE_DEPRECATED_OS=1 
export HDF5_USE_FILE_LOCKING=FALSE 


{} vdj --fastqs {} --id {} --sample {} --reference {} --jobmode /scratch/groups/hheyn/software/cellranger/6.1.1/external/martian/jobmanagers/slurm.template;
""".format( gem_id, cfg.cellranger_4_path, fastq_path, gem_id, gem_id, cfg.vdj_reference_path)
    job_script_file.write(job_script)
    job_script_file.close()
