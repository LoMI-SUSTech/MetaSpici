rule getSingle_genome:
    input:
        "assembly/{sample}_final.contigs.fa"
    output:
        directory("assembly/{sample}_single_genome/split_contigs")
    log:
        "assembly/{sample}_single_genome.log"
    params:
        threads=config.get("threads", 1),
        min_length=1000000,
        outdir="assembly/{sample}_single_genome"
    conda: "MetaSpicii_version1"
    shell:
        """
        mkdir -p {params.outdir}
        seqkit seq {input} --min-len {params.min_length} -j {params.threads} -o {params.outdir}/all_large_contigs.fa 2>> {log}
        echo "Extracted contigs > 1M:" >> {log}
        
        cd {params.outdir}
        seqkit split all_large_contigs.fa \
            --by-id \
            -j {params.threads} \
            -O split_contigs 2>> ../../{log}
        """

