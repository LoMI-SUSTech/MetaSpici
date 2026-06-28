rule LongRead_assembly:
    input:
        long_reads = config["long_reads"]["file"]
    output:
        "assembly/{sample}_final.contigs.fa"
    log:
        "assembly/{sample}_assembly.log"
    params:
        outdir = "assembly/{sample}_metaflye",
        threads = config.get("threads", 1)
    conda: "MetaSpicii_version1"
    shell:
        """
        long_type="{config[long_reads][type]}"
        
        if [[ $long_type == "pacbio" ]]; then
            flye --pacbio-raw {input.long_reads} \
                --out-dir {params.outdir} \
                --threads {params.threads} --meta >> {log} 2>&1
            mv {params.outdir}/assembly.fasta {output}
        
        elif [[ $long_type == "nanopore" ]]; then
            flye --nano-raw {input.long_reads} \
                --out-dir {params.outdir} \
                --threads {params.threads} --meta >> {log} 2>&1
            mv {params.outdir}/assembly.fasta {output}
                
        elif [[ $long_type == "pacbioHiFi" ]]; then
            metaMDBG asm \
                --out-dir {params.outdir} \
                --in-hifi {input.long_reads} \
                --threads {params.threads} >> {log} 2>&1
            gunzip -c {params.outdir}/contigs.fasta.gz > {output}
        
        elif [[ $long_type == "nanoporeR10" ]]; then
            metaMDBG asm \
                --out-dir {params.outdir} \
                --in-ont {input.long_reads} \
                --threads {params.threads} >> {log} 2>&1
            gunzip -c {params.outdir}/contigs.fasta.gz > {output}

        else
            echo "Unknown long read type: $long_type"
            exit 1
        fi
        """
