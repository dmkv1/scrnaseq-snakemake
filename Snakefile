import pandas as pd
import os


configfile: "config.yaml"


# Load samples table
samples = pd.read_excel("input.xlsx")
samples = samples.fillna("")

# Make samples and config available to the common module
import workflow.scripts.common as common

common.samples = samples
common.config = config

# Import helper functions
from workflow.scripts.common import *


# Include rules
include: "workflow/rules/cellranger.smk"


# Define samples with available data
gex_samples = samples["Sample"].tolist()
bcr_samples = samples[(samples["BCR_fq1"] != "") & (samples["BCR_fq2"] != "")][
    "Sample"
].tolist()
tcr_samples = samples[(samples["TCR_fq1"] != "") & (samples["TCR_fq2"] != "")][
    "Sample"
].tolist()


# Wildcard constraints
wildcard_constraints:
    sample="[^/]+",  # Match anything except slashes


rule all:
    input:
        # Gene expression outputs
        [
            f"results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5"
            for sample in gex_samples
        ],
