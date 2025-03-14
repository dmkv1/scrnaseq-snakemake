rule create_sce:
    input:
        stamp="stamps/cellranger/{sample}.stamp",
    output:
        sce="results/{sample}/sce/{sample}_sce_1_raw.rds",
    params:
        sample="{sample}",
        threads=config["resources"]["cellranger"]["threads"],
        counts_dir="results/{sample}/cellranger/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix",
        vdj_b="results/{sample}/cellranger/outs/per_sample_outs/{sample}/vdj_b/filtered_contig_annotations.csv",
        vdj_t="results/{sample}/cellranger/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv",
    conda:
        "../envs/scrnaseq_analysis.yaml"
    log:
        "logs/{sample}_R_create_sce.log",
    script:
        "../scripts/create_sce.R"


rule sce_QC:
    input:
        sce="results/{sample}/sce/{sample}_sce_1_raw.rds",
    output:
        sce_unfiltered="results/{sample}/sce/{sample}_sce_2_QC_unfiltered.rds",
        sce_filtered="results/{sample}/sce/{sample}_sce_2_QC_filtered.rds",
        plot_transcripts="results/{sample}/plots/{sample}_QC_transcripts.png",
        plot_doublets="results/{sample}/plots/{sample}_QC_doublets.png",
        plot_UMAPs="results/{sample}/plots/{sample}_QC_UMAPs.png",
    params:
        sample="{sample}",
        threads=config["resources"]["cellranger"]["threads"],
        random_seed=config["params"]["random_seed"],
        KNN=config["params"]["QC"]["KNN"],
        mito_perc=config["params"]["QC"]["mito_perc"],
        cluster_perc=config["params"]["QC"]["cluster_perc"],
    conda:
        "../envs/scrnaseq_analysis.yaml"
    log:
        "logs/{sample}_R_sce_QC.log",
    script:
        "../scripts/sce_QC.R"


rule sce_annotations:
    input:
        sce="results/{sample}/sce/{sample}_sce_2_QC_filtered.rds",
    output:
        sce_annotated="results/{sample}/sce/{sample}_sce_3_annotated.rds",
        plot_UMAPs="results/{sample}/plots/{sample}_QC_UMAPs.png",
        markers="results/{sample}/plots/{sample}_markers.rds",
    params:
        sample="{sample}",
        threads=config["resources"]["cellranger"]["threads"],
        random_seed=config["params"]["random_seed"],
        KNN=config["params"]["annotation"]["KNN"],
    conda:
        "../envs/scrnaseq_analysis.yaml"
    log:
        "logs/{sample}_R_sce_annotation.log",
    script:
        "../scripts/sce_annotation.R"
