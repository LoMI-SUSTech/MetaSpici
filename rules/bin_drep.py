import os
def has_short_reads(sample):
    return os.path.exists(f"qc/{sample}_clean_1.fastq")

if has_short_reads(config["sample"]):
    rule bin_drep:
        input:
            hic="binning/{sample}_final_hicbinset",
            large_contigs=lambda w: f"assembly/{w.sample}_single_genome/split_contigs" if config.get("long_reads") else [],
            short="binning/{sample}_final_shortbinset" 
        output:
            directory("binning/{sample}_final_binset")
        log:
            "binning/{sample}_drep.log"
        params:
            threads=lambda w: config.get("threads", 1),
            outdir="binning/{sample}_drep"
        conda: "MetaSpicii_version1"
        shell:
            """
            mkdir -p binning/{wildcards.sample}_combine_bins

            find {input.hic} -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                newname=$(basename "$file" .fa)
                newname=$(basename "$newname" .fasta)
                cp "$file" "binning/{wildcards.sample}_combine_bins/hic_$newname.fa"
            done
            
            if [ -n "{input.large_contigs}" ] && [ -d "{input.large_contigs}" ]; then
                find {input.large_contigs} -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                    newname=$(basename "$file" .fa)
                    newname=$(basename "$newname" .fasta)
                    cp "$file" "binning/{wildcards.sample}_combine_bins/large_$newname.fa"
                done
            fi

            find "{input.short}" -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                newname=$(basename "$file" .fa)
                newname=$(basename "$newname" .fasta)
                cp "$file" "binning/{wildcards.sample}_combine_bins/short_$newname.fa"
            done

            dRep dereplicate {params.outdir} \
                -g binning/{wildcards.sample}_combine_bins/*.fa \
                -comp 50 -con 5 -sa 0.99 \
                -p {params.threads} >> {log} 2>&1

            mv {params.outdir}/dereplicated_genomes {output}
            echo "dRep with short reads completed" >> {log}
            """
            
else:
    rule bin_drep:
        input:
            hic="binning/{sample}_final_hicbinset",
            large_contigs=lambda w: f"assembly/{w.sample}_single_genome/split_contigs" if config.get("long_reads") else []
        output:
            directory("binning/{sample}_final_binset")
        log:
            "binning/{sample}_drep.log"
        params:
            threads=lambda w: config.get("threads", 1),
            outdir="binning/{sample}_drep"
        conda: "MetaSpicii_version1"
        shell:
            """
            mkdir -p binning/{wildcards.sample}_combine_bins

            find {input.hic} -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                newname=$(basename "$file" .fa)
                newname=$(basename "$newname" .fasta)
                cp "$file" "binning/{wildcards.sample}_combine_bins/hic_$newname.fa"
            done
            
            if [ -n "{input.large_contigs}" ] && [ -d "{input.large_contigs}" ]; then
                find {input.large_contigs} -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                    newname=$(basename "$file" .fa)
                    newname=$(basename "$newname" .fasta)
                    cp "$file" "binning/{wildcards.sample}_combine_bins/large_$newname.fa"
                done
            fi

            dRep dereplicate {params.outdir} \
                -g binning/{wildcards.sample}_combine_bins/*.fa \
                -comp 50 -con 5 -sa 0.99 \
                -p {params.threads} >> {log} 2>&1

            mv {params.outdir}/dereplicated_genomes {output}
            echo "dRep without short reads completed" >> {log}
            """
