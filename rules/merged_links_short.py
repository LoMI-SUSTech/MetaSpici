rule merge_MAGs:
    input:
        "binning/{sample}_final_binset"
    output:
        "mapping/reference/{sample}_MAG_merged.fasta"
    shell:
        """
        mkdir -p mapping/reference

        for f in {input}/*.fa; do
            name=$(basename $f .fa)
            awk -v name="$name" '
                BEGIN {{seq=""}}
                /^>/ {{next}}
                {{seq=seq$0}}
                END {{print ">"name"\\n"seq}}
            ' $f
        done > {output}
        """

rule prepare_MGE:
    input:
        "output/MGE_detect/metator_mge_{sample}/mges_binned.fa"
    output:
        "mapping/reference/{sample}_MGE.fasta"
    shell:
        """
        mkdir -p mapping/reference
        cp {input} {output}
        """

rule build_reference:
    input:
        mag="mapping/reference/{sample}_MAG_merged.fasta",
        mge="mapping/reference/{sample}_MGE.fasta"
    output:
        "mapping/reference/{sample}_combined.fasta"
    shell:
        """
        cat {input.mag} {input.mge} > {output}
        """


rule metaxa_rrna:
    input:
        "mapping/reference/{sample}_combined.fasta"
    output:
        ssu="mapping/reference/reference_masked/{sample}.metaxa2_ssu.extraction.results",
        lsu="mapping/reference/reference_masked/{sample}.metaxa2_lsu.extraction.results",
        bed="mapping/reference/reference_masked/{sample}_masked_metaxa.bed"
    threads: config["threads"]
    conda:
        "metaxa-2.2"
    shell:
        r"""
        set -euo pipefail

        mkdir -p mapping/reference/reference_masked

        metaxa2 --plus --mode m --cpu {threads} --multi_thread T --table T -g ssu --not_found T \
            -i {input} -o mapping/reference/reference_masked/{wildcards.sample}.metaxa2_ssu

        metaxa2 --plus --mode m --cpu {threads} --multi_thread T --table T -g lsu --not_found T \
            -i {input} -o mapping/reference/reference_masked/{wildcards.sample}.metaxa2_lsu

        cut -f 1,9,10 mapping/reference/reference_masked/{wildcards.sample}.metaxa2_ssu.extraction.results > {output.bed}
        cut -f 1,9,10 mapping/reference/reference_masked/{wildcards.sample}.metaxa2_lsu.extraction.results >> {output.bed}
        """

rule barrnap_rrna:
    input:
        "mapping/reference/{sample}_combined.fasta"
    output:
        "mapping/reference/reference_masked/{sample}_masked_barrnap.bed"
    threads: config["threads"]
    conda:
        "barnap"
    shell:
        """
        barrnap --kingdom bac --threads {threads} --reject 0.3 {input} > tmp_bac.gff
        barrnap --kingdom arc --threads {threads} --reject 0.3 {input} > tmp_arc.gff
        barrnap --kingdom euk --threads {threads} --reject 0.3 {input} > tmp_euk.gff

        cut -f 1,4,5 tmp_bac.gff > {output}
        cut -f 1,4,5 tmp_arc.gff >> {output}
        cut -f 1,4,5 tmp_euk.gff >> {output}

        rm tmp_bac.gff tmp_arc.gff tmp_euk.gff
        """

rule trna_mask:
    input:
        "mapping/reference/{sample}_combined.fasta"
    output:
        "mapping/reference/reference_masked/{sample}_masked_trna.bed"
    threads: config["threads"]
    conda:
        "trnascan-se"
    shell:
        """
        tRNAscan-SE -A -b tmp_archaea.bed --thread {threads} {input}
        tRNAscan-SE -B -b tmp_bacteria.bed --thread {threads} {input}

        cat tmp_archaea.bed tmp_bacteria.bed | cut -f 1,2,3 > {output}

        rm tmp_archaea.bed tmp_bacteria.bed
        """

rule dust_mask:
    input:
        "mapping/reference/{sample}_combined.fasta"
    output:
        "mapping/reference/reference_masked/{sample}_masked_dust.bed"
    conda:
        "blast"
    shell:
        """
        dustmasker -in {input} -out tmp.lowcom.out -outfmt acclist
        cut -d '>' -f 2 tmp.lowcom.out > {output}
        rm tmp.lowcom.out
        """

rule apply_mask:
    input:
        fasta="mapping/reference/{sample}_combined.fasta",
        metaxa="mapping/reference/reference_masked/{sample}_masked_metaxa.bed",
        barrnap="mapping/reference/reference_masked/{sample}_masked_barrnap.bed",
        trna="mapping/reference/reference_masked/{sample}_masked_trna.bed",
        dust="mapping/reference/reference_masked/{sample}_masked_dust.bed"
    output:
        "mapping/reference/reference_masked/{sample}_masked.fasta"
    conda:
        "map"
    shell:
        """
        cat {input.metaxa} {input.barrnap} {input.trna} {input.dust} > tmp_mask.bed

        awk -F'\t' '$2>=0' tmp_mask.bed > tmp_mask_clean.bed

        bedtools maskfasta \
            -fi {input.fasta} \
            -bed tmp_mask_clean.bed \
            -fo {output}

        rm tmp_mask.bed tmp_mask_clean.bed
        """
