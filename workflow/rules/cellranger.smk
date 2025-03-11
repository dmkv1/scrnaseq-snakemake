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
        time=config["resources"]["cellranger"]["time"],
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
            > {log} 2>&1
        
        mv {params.sample}/outs/* results/cellranger/{params.sample}/outs/
        rm -rf {params.sample}
        """


# Rule for Cell Ranger VDJ (BCR analysis)
rule cellranger_vdj_bcr:
    input:
        r1=lambda wildcards: get_bcr_fastqs(wildcards)["r1"],
        r2=lambda wildcards: get_bcr_fastqs(wildcards)["r2"],
    output:
        annotations="results/cellranger/{sample}/outs/vdj_b/filtered_contig_annotations.csv",
        contigs="results/cellranger/{sample}/outs/vdj_b/filtered_contig.fasta",
    params:
        cellranger=config["tools"]["cellranger"],
        sample="{sample}",
        fastq_dir=lambda wildcards: get_fastq_dir(wildcards, "BCR"),
        reference=lambda _: config["reference"]["vdj"],
    threads: config["resources"]["cellranger"]["threads"]
    resources:
        mem_mb=config["resources"]["cellranger"]["memory"],
        time=config["resources"]["cellranger"]["time"],
    log:
        "logs/cellranger/{sample}_vdj_bcr.log",
    shell:
        """
        mkdir -p results/cellranger/{params.sample}/outs/vdj_b
        
        {params.cellranger} vdj \
            --id={params.sample}_bcr \
            --reference={params.reference} \
            --fastqs={params.fastq_dir} \
            --sample={params.sample} \
            --localcores={threads} \
            --localmem={resources.mem_mb} \
            > {log} 2>&1
        
        mv {params.sample}_bcr/outs/* results/cellranger/{params.sample}/outs/vdj_b/
        rm -rf {params.sample}_bcr
        """


# Rule for Cell Ranger VDJ (TCR analysis)
rule cellranger_vdj_tcr:
    input:
        r1=lambda wildcards: get_tcr_fastqs(wildcards)["r1"],
        r2=lambda wildcards: get_tcr_fastqs(wildcards)["r2"],
    output:
        annotations="results/cellranger/{sample}/outs/vdj_t/filtered_contig_annotations.csv",
        contigs="results/cellranger/{sample}/outs/vdj_t/filtered_contig.fasta",
    params:
        cellranger=config["tools"]["cellranger"],
        sample="{sample}",
        fastq_dir=lambda wildcards: get_fastq_dir(wildcards, "TCR"),
        reference=lambda _: config["reference"]["vdj"],
    threads: config["resources"]["cellranger"]["threads"]
    resources:
        mem_mb=config["resources"]["cellranger"]["memory"],
        time=config["resources"]["cellranger"]["time"],
    log:
        "logs/cellranger/{sample}_vdj_tcr.log",
    shell:
        """
        mkdir -p results/cellranger/{params.sample}/outs/vdj_t
        
        {params.cellranger} vdj \
            --id={params.sample}_tcr \
            --reference={params.reference} \
            --fastqs={params.fastq_dir} \
            --sample={params.sample} \
            --localcores={threads} \
            --localmem={resources.mem_mb} \
            > {log} 2>&1
        
        mv {params.sample}_tcr/outs/* results/cellranger/{params.sample}/outs/vdj_t/
        rm -rf {params.sample}_tcr
        """
