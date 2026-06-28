rule HiC_binning:
    input:
        "assembly/{sample}_final_min2k.contigs.fa",
        "qc/{sample}_hic_clean_1.fastq",
        "qc/{sample}_hic_clean_2.fastq"
    output:
        directory("binning/{sample}_final_hicbinset")
    log:
        "binning/{sample}_hicbinning.log"
    params:
        threads=config.get("threads", 1),
        outdir="binning/{sample}_metator",
        enzyme=config.get("enzyme")
    conda: "MetaSpicii_version1"
    shell:
        """
        export CHECKVDB=/lomi_home/zichun/MetaSpicii/db/checkv-db-v1.5/
        metator pipeline -1 {input[1]} \
            -2 {input[2]} \
            -a {input[0]} \
            -e {params.enzyme} \
            -o {params.outdir} \
            -t {params.threads} \
            >> {log} 2>&1
        mv {params.outdir}/final_bin_unscaffold {output}
        """
