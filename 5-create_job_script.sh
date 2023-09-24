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
