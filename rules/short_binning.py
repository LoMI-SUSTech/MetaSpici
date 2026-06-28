rule short_binning:
    input:
        contigs="assembly/{sample}_final_min2k.contigs.fa",
        r1="qc/{sample}_clean_1.fastq", 
        r2="qc/{sample}_clean_2.fastq"
    output:
        directory("binning/{sample}_final_shortbinset")
    log:
        "binning/{sample}_basalt/{sample}_binning.log"
    params:
        threads=lambda wildcards: config.get("threads", 1),
        outdir="binning/{sample}_basalt"
    conda: "BASALT"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        ln -sf ../../{input.contigs} .
        ln -sf ../../{input.r1} .
        ln -sf ../../{input.r2} .

        contig_basename=$(basename {input.contigs})
        r1_basename=$(basename {input.r1})
        r2_basename=$(basename {input.r2})
        
        BASALT -a "$contig_basename" -t {params.threads} \
               -s "$r1_basename,$r2_basename" \
               --mode continue \
               --module autobinning >> {wildcards.sample}_binning.log 2>&1
        
        BASALT -a "$contig_basename" -t {params.threads} \
               -s "$r1_basename,$r2_basename" \
               --mode continue \
               --module refinement >> {wildcards.sample}_binning.log 2>&1
        
        mv BestBinset_outlier_refined_filtrated_retrieved_retrieved ../{wildcards.sample}_final_shortbinset
        """
