rule length_filtering:
    input:
        contigs="assembly/{sample}_final.contigs.fa"
    output:
        output_2k="assembly/{sample}_final_min2k.contigs.fa"
    conda: "MetaSpicii_version1"
    shell:
        """
        seqmagick convert --min-length 2000 {input.contigs} {output.output_2k}
        """
