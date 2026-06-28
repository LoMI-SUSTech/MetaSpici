rule hic_cutsite:
    input:
        r1="qc/{sample}_hic_clean_1.fastq",
        r2="qc/{sample}_hic_clean_2.fastq"
    output:
        r1="qc/{sample}_hiccutsite_clean_1.fastq",
        r2="qc/{sample}_hiccutsite_clean_2.fastq"
    log:
        "qc/{sample}_hic_cutsite.log"
    params:
        threads=config.get("threads", 1),
        enzyme=config.get("enzyme"),
        out_prefix="{sample}_hiccutsite"
    conda: "MetaSpicii_version1"
    shell:
        """
        hicstuff cutsite -1 {input.r1} -2 {input.r2} -p {params.out_prefix} \
            -e {params.enzyme} -s 30 -t {params.threads} >> {log} 2>&1
        gunzip -c {params.out_prefix}_R1.fq.gz > {output.r1}  && rm {params.out_prefix}_R1.fq.gz
        gunzip -c {params.out_prefix}_R2.fq.gz > {output.r2}  && rm {params.out_prefix}_R2.fq.gz
        """


if "short_reads" not in config:
    print("INFO: No short reads provided, creating symlinks from HiC reads for assembly/binning steps.")
    sample_name = config["sample"]
    rule create_assembly_input_symlinks:
        input:
            r1=f"qc/{sample_name}_hiccutsite_clean_1.fastq",
            r2=f"qc/{sample_name}_hiccutsite_clean_2.fastq"
        output:
            r1=f"qc/{sample_name}_clean_1.fastq",
            r2=f"qc/{sample_name}_clean_2.fastq"
        shell:
            """
            ln -sf $(realpath {input.r1}) {output.r1}
            ln -sf $(realpath {input.r2}) {output.r2}
            """
