# Cell Ranger rules for scRNA-seq and VDJ analysis


# Rule for Cell Ranger count (GEX analysis)
rule cellranger_count:
    input:
        r1=lambda wildcards: get_gex_fastqs(wildcards)["r1"],
        r2=lambda wildcards: get_gex_fastqs(wildcards)["r2"],
    output:
        matrix="results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
        summary="results/cellranger/{sample}/outs/metrics_summary.csv",
    params:
        cellranger=config["tools"]["cellranger"],
        sample="{sample}",
        fastq_dir=lambda wildcards: get_fastq_dir(wildcards, "GEX"),
        expected_cells=get_expected_cells,
        transcriptome=lambda _: config["reference"]["transcriptome"],
    threads: config["resources"]["cellranger"]["threads"]
    resources:
        mem_mb=config["resources"]["cellranger"]["memory"],
    log:
        "logs/cellranger/{sample}_count.log",
    shell:
        """
        mkdir -p results/cellranger/{params.sample}/outs
        
        {params.cellranger} count \
            --id={params.sample} \
            --transcriptome={params.transcriptome} \
            --fastqs={params.fastq_dir} \
            --sample={params.sample} \
            --expect-cells={params.expected_cells} \
            --localcores={threads} \
            --localmem={resources.mem_mb} \
            --create-bam true \
            > {log} 2>&1
        
        mv {params.sample}/outs/* results/cellranger/{params.sample}/outs/
        rm -rf {params.sample}
        """
