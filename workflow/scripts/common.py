import os

# These variables will be populated by the main Snakefile
samples = None
config = None

def get_gex_fastqs(wildcards):
    """Return the fastq paths for gene expression data"""
    if samples is None:
        raise ValueError("samples DataFrame is not initialized")
    
    matching_rows = samples[samples["Sample"] == wildcards.sample]
    if len(matching_rows) == 0:
        raise ValueError(f"No sample found with name '{wildcards.sample}'. Available samples: {samples['Sample'].tolist()}")
    
    sample_row = matching_rows.iloc[0]
    return {
        "r1": sample_row["GEX_fq1"],
        "r2": sample_row["GEX_fq2"]
    }

def get_bcr_fastqs(wildcards):
    """Return the fastq paths for BCR data"""
    sample_row = samples[samples["Sample"] == wildcards.sample].iloc[0]
    return {
        "r1": sample_row["BCR_fq1"],
        "r2": sample_row["BCR_fq2"]
    }

def get_tcr_fastqs(wildcards):
    """Return the fastq paths for TCR data"""
    sample_row = samples[samples["Sample"] == wildcards.sample].iloc[0]
    return {
        "r1": sample_row["TCR_fq1"],
        "r2": sample_row["TCR_fq2"]
    }

def get_fastq_dir(wildcards, data_type="GEX"):
    """
    Return the directory containing the FASTQ files for a given sample and data type.
    """
    sample_row = samples[samples["Sample"] == wildcards.sample].iloc[0]
    return os.path.dirname(sample_row[f"{data_type}_fq1"])

def get_expected_cells(wildcards):
    """
    Return the expected number of cells for a given sample.
    """
    return samples[samples["Sample"] == wildcards.sample]["expected_cells"].iloc[0]

def get_sample_paths(wildcards):
    """
    Return all file paths for a given sample.
    """
    sample_row = samples[samples["Sample"] == wildcards.sample].iloc[0]
    
    result = {}
    
    # Add GEX paths if they exist
    if sample_row["GEX_fq1"] and sample_row["GEX_fq2"]:
        result["gex_r1"] = sample_row["GEX_fq1"]
        result["gex_r2"] = sample_row["GEX_fq2"]
    
    # Add BCR paths if they exist
    if sample_row["BCR_fq1"] and sample_row["BCR_fq2"]:
        result["bcr_r1"] = sample_row["BCR_fq1"]
        result["bcr_r2"] = sample_row["BCR_fq2"]
    
    # Add TCR paths if they exist
    if sample_row["TCR_fq1"] and sample_row["TCR_fq2"]:
        result["tcr_r1"] = sample_row["TCR_fq1"]
        result["tcr_r2"] = sample_row["TCR_fq2"]
    
    return result