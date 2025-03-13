import os
import re


def extract_fastq_prefix(filepath):
    filename = os.path.basename(filepath)
    
    # 10x Genomics naming format: PREFIX_S1_L001_R1_001.fastq.gz
    match = re.match(r'^(.+?)_S\d+_L\d{3}_R\d_\d{3}\.fastq\.gz$', filename)
    if match:
        return match.group(1)
    
    # If regular expression fails, raise an error
    raise ValueError(f"Fastq file '{filename}' does not match the expected 10X Genomics format: PREFIX_S1_L001_R1_001.fastq.gz")


def create_multi_config(sample, samples_df, config, output_file):
    sample_row = samples_df[samples_df["Sample"] == sample].iloc[0]

    has_gex = sample_row["GEX_fq1"] != "" and sample_row["GEX_fq2"] != ""
    has_bcr = sample_row["BCR_fq1"] != "" and sample_row["BCR_fq2"] != ""
    has_tcr = sample_row["TCR_fq1"] != "" and sample_row["TCR_fq2"] != ""

    csv_lines = []

    if has_gex:
        csv_lines.extend(
            [
                "[gene-expression]",
                f"reference,{config['reference']['transcriptome']}",
                f"expect-cells,{sample_row['expected_cells']}",
                "create-bam,true",
                "",
            ]
        )

    if has_bcr or has_tcr:
        csv_lines.extend(["[vdj]", f"reference,{config['reference']['vdj']}", ""])

    csv_lines.append("[libraries]")
    csv_lines.append("fastq_id,fastqs,feature_types")

    if has_gex:
        gex_dir = os.path.dirname(sample_row["GEX_fq1"])
        gex_prefix = extract_fastq_prefix(sample_row["GEX_fq1"])
        csv_lines.append(f"{gex_prefix},{gex_dir},gene expression")

    if has_bcr:
        bcr_dir = os.path.dirname(sample_row["BCR_fq1"])
        bcr_prefix = extract_fastq_prefix(sample_row["BCR_fq1"])
        csv_lines.append(f"{bcr_prefix},{bcr_dir},VDJ-B")

    if has_tcr:
        tcr_dir = os.path.dirname(sample_row["TCR_fq1"])
        tcr_prefix = extract_fastq_prefix(sample_row["TCR_fq1"])
        csv_lines.append(f"{tcr_prefix},{tcr_dir},VDJ-T")

    with open(output_file, "w") as f:
        f.write("\n".join(csv_lines))