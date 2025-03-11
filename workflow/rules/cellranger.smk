rule cellranger_count:
    input:
        r1=lambda wildcards: get_gex_fastqs(wildcards)["r1"],
        r2=lambda wildcards: get_gex_fastqs(wildcards)["r2"],
    output:
        matrix="results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
        summary="results/cellranger/{sample}/outs/metrics_summary.csv",
    params:
        cellranger=config["tools"]["cellranger"],
        transcriptome=lambda _: config["reference"]["transcriptome"],
        outdir="results/cellranger/{sample}",
        sample="{sample}",
        fastq_dir=lambda wildcards: get_fastq_dir(wildcards, "GEX"),
        sample_prefix=lambda wildcards: get_fastq_prefix(wildcards, "GEX"),
        expected_cells=get_expected_cells,
    threads: config["resources"]["cellranger"]["threads"]
    resources:
        mem_gb=config["resources"]["cellranger"]["memory"],
    log:
        "logs/cellranger/{sample}_count.log",
    shell:
        """
        mkdir -p {params.outdir}
        
        {params.cellranger} count \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            --disable-ui \
            --create-bam true \
            --nosecondary true \
            --transcriptome={params.transcriptome} \
            --output-dir {params.outdir} \
            --id={params.sample} \
            --fastqs={params.fastq_dir} \
            --sample={params.sample_prefix} \
            --expect-cells={params.expected_cells} \
            > {log} 2>&1
        """
