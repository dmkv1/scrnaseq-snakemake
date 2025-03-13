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
        directory("results/{sample}/cellranger/outs"),
        counts="results/{sample}/cellranger/outs/per_sample_outs/{sample}/sample_filtered_feature_bc_matrix/matrix.mtx.gz",
        features="results/{sample}/cellranger/outs/per_sample_outs/{sample}/sample_filtered_feature_bc_matrix/features.tsv.gz",
        barcodes="results/{sample}/cellranger/outs/per_sample_outs/{sample}/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
        vdj_b_csv="results/{sample}/cellranger/outs/per_sample_outs/{sample}/vdj_b/filtered_contig_annotations.csv",
        vdj_t_csv="results/{sample}/cellranger/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv",
        websummary="results/{sample}/cellranger/outs/per_sample_outs/{sample}/web_summary.html",
        counter=temp(directory("results/{sample}/cellranger/SC_MULTI_CS")),
    params:
        cellranger=config["tools"]["cellranger"],
        sample="{sample}",
        rundir="results/{sample}/cellranger_run",
        outdir="results/{sample}/cellranger",
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

        mv {params.rundir} {params.outdir}
        """
