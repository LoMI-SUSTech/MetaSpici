rule bin_drep:
    input:
        hic="binning/{sample}_final_hicbinset",
        short="binning/{sample}_final_shortbinset"
    output:
        directory("binning/{sample}_final_binset")
    log:
        "binning/{sample}_drep.log"
    params:
        threads=lambda w: config.get("threads", 1),
        outdir="binning/{sample}_drep",
        mode=config["mode"],
        has_long=str(config.get("long_reads", False)).lower()
    conda:
        "MetaSpicii_version1"
    shell:
        r"""
        mkdir -p binning/{wildcards.sample}_combine_bins

        echo "MODE: {params.mode}" >> {log}

        # --------------------------
        # 1️⃣ HIC（永远存在）
        # --------------------------
        if [ -d "{input.hic}" ]; then
            find {input.hic} -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                newname=$(basename "$file" .fa)
                newname=$(basename "$newname" .fasta)
                cp "$file" "binning/{wildcards.sample}_combine_bins/hic_$newname.fa"
            done
        fi

        # --------------------------
        # 2️⃣ short（只有非 long_hic 才用）
        # --------------------------
        if [ "{params.mode}" != "long_hic" ]; then
            if [ -d "{input.short}" ]; then
                find {input.short} -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                    newname=$(basename "$file" .fa)
                    newname=$(basename "$newname" .fasta)
                    cp "$file" "binning/{wildcards.sample}_combine_bins/short_$newname.fa"
                done
            fi
        fi

        # --------------------------
        # 3️⃣ long（可选）
        # --------------------------
        if [ "{params.has_long}" = "true" ]; then
            long_dir="assembly/{wildcards.sample}_single_genome/split_contigs"
            if [ -d "$long_dir" ]; then
                find "$long_dir" -name "*.fa" -o -name "*.fasta" 2>/dev/null | while read file; do
                    newname=$(basename "$file" .fa)
                    newname=$(basename "$newname" .fasta)
                    cp "$file" "binning/{wildcards.sample}_combine_bins/long_$newname.fa"
                done
            fi
        fi

        # --------------------------
        # 4️⃣ dRep
        # --------------------------
        dRep dereplicate {params.outdir} \
            -g binning/{wildcards.sample}_combine_bins/*.fa \
            -comp 50 -con 5 -sa 0.99 \
            -p {params.threads} >> {log} 2>&1

        mv {params.outdir}/dereplicated_genomes {output}

        echo "dRep finished (mode={params.mode})" >> {log}
        """
