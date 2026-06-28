rule hic_read_qc:
    input:
        r1 = lambda wildcards: config["hic_reads"]["r1"],
        r2 = lambda wildcards: config["hic_reads"]["r2"]
    output:
        f"qc/{config['sample']}_hic_clean_1.fastq",
        f"qc/{config['sample']}_hic_clean_2.fastq",
        f"qc/{config['sample']}_hic.json",
        f"qc/{config['sample']}_hic.html"
    log:
        f"qc/{config['sample']}_hic_qc.log"
    params:
        threads=config.get("threads", 1)
    conda: "MetaSpicii_version1"
    shell:
        """
            fastp -i {input.r1} -I {input.r2} \
            -o {output[0]} -O {output[1]} \
            -w {params.threads} -D -q 30 -l 30 \
            --cut_front --cut_tail \
            -r -M 30 -j {output[2]} -h {output[3]} \
            --umi_delim / >> {log} 2>&1
        """
