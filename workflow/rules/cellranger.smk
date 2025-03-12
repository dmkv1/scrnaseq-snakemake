rule cellranger_count:
    input:
        r1=lambda wildcards: get_gex_fastqs(wildcards)["r1"],
        r2=lambda wildcards: get_gex_fastqs(wildcards)["r2"],
    output:
        bam="results/cellranger/{sample}/outs/possorted_genome_bam.bam",
        matrix.gz="results/cellranger/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        features.gz="results/cellranger/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes.gz="results/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        raw.matrix.gz="results/cellranger/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz",
        raw.features.gz="results/cellranger/{sample}/outs/raw_feature_bc_matrix/features.tsv.gz",
        raw.barcodes.gz="results/cellranger/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
        summary="results/cellranger/{sample}/outs/web_summary.html",
        counter=temp(directory("results/cellranger/{sample}/SC_RNA_COUNTER_CS")),
    params:
        cellranger=config["tools"]["cellranger"],
        transcriptome=:config["reference"]["transcriptome"],
        outdir="results/cellranger/{sample}",
        sample="{sample}",
        fastq_dir=lambda wildcards: get_fastq_dir(wildcards, "GEX"),
        sample_prefix=lambda wildcards: get_fastq_prefix(wildcards, "GEX"),
        expected_cells=get_expected_cells,
    threads: config["resources"]["cellranger"]["threads"],
    resources:
        mem_gb=config["resources"]["cellranger"]["memory"],
    log:
        "logs/cellranger/{sample}_count.log",
    shell:
        """
        rm -rf {params.outdir}

        {params.cellranger} count \
            --id={params.sample} \
            --create-bam true \
            --transcriptome={params.transcriptome} \
            --fastqs={params.fastq_dir} \
            --sample={params.sample_prefix} \
            --expect-cells={params.expected_cells} \
            --output-dir {params.outdir} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            --disable-ui \
            --nosecondary \
            > {log} 2>&1
        """

