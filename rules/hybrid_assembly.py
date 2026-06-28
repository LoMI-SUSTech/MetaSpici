#long reads assembly with metaflye
rule Hybrid_assembly:
    input:
        short_reads_r1 =  "qc/{sample}_clean_1.fastq",
        short_reads_r2 = "qc/{sample}_clean_2.fastq",
        Long_reads = config["long_reads"]["file"]
    output:
        "assembly/{sample}_hybrid_contigs.fa"
    log:
        "assembly/{sample}_assembly.log"
    params:
        outdir="assembly/{sample}_operams",
        threads=config.get("threads", 1)
    conda: "MetaSpicii_version1"
    shell:
        """
        perl /lomi_home/sunzichun/metaspicii_realsamples/scripts/OPERA-MS.pl \
            --short-read1 {input[0]} \
            --short-read2 {input[1]} \
            --long-read {input[2]} \
            --no-ref-clustering \
            --num-processors {params.threads} \
            --out-dir {params.outdir} >> {log} 2>&1
        mv {params.outdir}/contigs.polished.fasta {output}
        """



rule LongRead_assembly:
    input:
        long_reads = config["long_reads"]["file"]
    output:
        "assembly/{sample}_long.contigs.fa"
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

rule merge_assemblies:
    input:
        hybrid_contigs = "assembly/{sample}_hybrid_contigs.fa",
        long_contigs = "assembly/{sample}_long.contigs.fa"
    output:
        "assembly/{sample}_final.contigs.fa"
    log:
        "assembly/{sample}_merge.log"
    params:
        outdir = "assembly/{sample}_quickmerge",
        threads = config.get("threads", 1)
    conda: "MetaSpicii_version1"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        ln -sf ../{wildcards.sample}_hybrid_contigs.fa {wildcards.sample}_hybrid_contigs.fa
        ln -sf ../{wildcards.sample}_long.contigs.fa {wildcards.sample}_long.contigs.fa 
        export PATH=/lomi_home/wenxiu/software/quickmerge:/lomi_home/wenxiu/software/quickmerge/MUMmer3.23:$PATH
        merge_wrapper.py -l 1000 {wildcards.sample}_hybrid_contigs.fa {wildcards.sample}_long.contigs.fa >> {wildcards.sample}_merge.log 2>&1
        mv merged_out.fasta ../{wildcards.sample}_final.contigs.fa
        """
