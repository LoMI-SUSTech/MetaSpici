rule short_assembly:
    input:
        "qc/{sample}_clean_1.fastq",
        "qc/{sample}_clean_2.fastq"
    output:
        "assembly/{sample}_final.contigs.fa"
    log:
        "assembly/{sample}_assembly.log"
    params:
        threads=config.get("threads", 1),
        memory=0.9,
        preset="meta-sensitive",
        outdir="assembly/{sample}_megahit"
    conda: "MetaSpicii_version1"
    shell:
        """
        rm -rf {params.outdir}
        megahit --presets {params.preset} \
            -t {params.threads} \
            -m {params.memory} \
            -1 {input[0]} \
            -2 {input[1]} \
            -o {params.outdir} >> {log} 2>&1
        mv {params.outdir}/final.contigs.fa {output}
        """
