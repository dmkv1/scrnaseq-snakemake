rule create_multi_config:
    output:
        config_csv="results/{sample}/multi_config.csv",
    params:
        sample="{sample}",
    run:
        from workflow.scripts.create_multi_config import create_multi_config

        create_multi_config(params.sample, samples_df, config, output.config_csv)


rule cellranger_multi:
    input:
        config_csv="results/{sample}/multi_config.csv",
    output:
        stamp="stamps/cellranger/{sample}.stamp",
    params:
        cellranger=config["tools"]["cellranger"],
        sample="{sample}",
        rundir="results/{sample}/cellranger",
    threads: config["resources"]["cellranger"]["threads"]
    resources:
        mem_gb=config["resources"]["cellranger"]["memory"],
    log:
        "logs/cellranger/{sample}_multi.log",
    shell:
        """
        {params.cellranger} multi \
            --id={params.sample} \
            --csv={input.config_csv} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            --output-dir {params.rundir} \
            --disable-ui \
            > {log} 2>&1
        touch {output.stamp}
        """


rule create_sce:
    input:
        stamp="stamps/cellranger/{sample}.stamp",
    output:
        sce="results/{sample}/sce/{sample}_sce_raw.rds",
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
