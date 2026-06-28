rule short_read_qc:
    input:
        r1 = lambda wildcards: config["short_reads"]["r1"],
        r2 = lambda wildcards: config["short_reads"]["r2"]
    output:
        f"qc/{config['sample']}_clean_1.fastq",
        f"qc/{config['sample']}_clean_2.fastq",
        f"qc/{config['sample']}.json",
        f"qc/{config['sample']}.html"
    log:
        f"qc/{config['sample']}_qc.log"
    params:
        threads=config.get("threads", 1)
    conda: "MetaSpicii_version1"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
            -o {output[0]} -O {output[1]} \
            -w {params.threads} -D -q 30 -l 75 \
            --cut_front --cut_tail \
            -r -M 30 -j {output[2]} -h {output[3]} \
            --umi_delim / >> {log} 2>&1
        """
