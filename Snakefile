import pandas as pd
import os


configfile: "config.yaml"


# Load samples table
samples_df = pd.read_excel("input.xlsx")
samples_df = samples_df.fillna("")


# Include rules
include: "workflow/rules/cellranger.smk"


# Define samples with available data
samples = samples_df["Sample"].tolist()


# Wildcard constraints
wildcard_constraints:
    sample="[^/]+",  # Match anything except slashes


rule all:
    input:
        [
            f"results/{sample}/cellranger/outs/per_sample_outs/{sample}/sample_filtered_feature_bc_matrix/matrix.mtx.gz"
            for sample in samples
        ],
