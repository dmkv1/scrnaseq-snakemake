import os

# These variables will be populated by the main Snakefile
samples = None
config = None


def get_gex_fastqs(wildcards):
    if samples is None:
        raise ValueError("samples DataFrame is not initialized")

    matching_rows = samples[samples["Sample"] == wildcards.sample]
    if len(matching_rows) == 0:
        raise ValueError(
            f"No sample found with name '{wildcards.sample}'. Available samples: {samples['Sample'].tolist()}"
        )

    sample_row = matching_rows.iloc[0]
    return {"r1": sample_row["GEX_fq1"], "r2": sample_row["GEX_fq2"]}


def get_fastq_dir(wildcards, data_type="GEX"):
    sample_row = samples[samples["Sample"] == wildcards.sample].iloc[0]
    return os.path.dirname(sample_row[f"{data_type}_fq1"])


def get_fastq_prefix(wildcards, data_type="GEX"):
    sample_row = samples[samples["Sample"] == wildcards.sample].iloc[0]
    fastq_path = sample_row[f"{data_type}_fq1"]
    filename = os.path.basename(fastq_path)
    parts = filename.split("_S")
    if len(parts) >= 1:
        return parts[0]
    return filename


def get_expected_cells(wildcards):
    return samples[samples["Sample"] == wildcards.sample]["expected_cells"].iloc[0]
