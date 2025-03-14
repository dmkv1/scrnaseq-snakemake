import pandas as pd
import os


configfile: "config.yaml"


# Load samples table
samples_df = pd.read_excel("input.xlsx")
samples_df = samples_df.fillna("")


# Include rules
include: "workflow/rules/cellranger.smk"
include: "workflow/rules/SCE_processing.smk"


# Define samples with available data
samples = samples_df["Sample"].tolist()


# Wildcard constraints
wildcard_constraints:
    sample="[^/]+",  # Match anything except slashes


rule all:
    input:
        [f"stamps/cellranger/{sample}.stamp" for sample in samples],
        [f"results/{sample}/sce/{sample}_sce_1_raw.rds" for sample in samples],
        [f"results/{sample}/sce/{sample}_sce_2_QC_filtered.rds" for sample in samples],
        [f"results/{sample}/sce/{sample}_sce_3_annotated.rds" for sample in samples],
