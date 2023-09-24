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
