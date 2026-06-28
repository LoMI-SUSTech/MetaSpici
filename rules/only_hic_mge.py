########################################
# 1️⃣ Contig length filtering (>5kb)
########################################

rule mge_length_filtering:
    input:
        contigs="assembly/{sample}_final.contigs.fa"
    output:
        "output/MGE_detect/{sample}.contigs.min5k.fa"
    conda:
        "MetaSpicii_version1"
    shell:
        """
        mkdir -p output/MGE_detect
        seqmagick convert --min-length 5000 {input.contigs} {output}
        """


########################################
# 2️⃣ geNomad MGE detection
########################################
rule genomad_detection:
    input:
        "output/MGE_detect/{sample}.contigs.min5k.fa"
    output:
        virus="output/MGE_detect/genomad_{sample}/{sample}.contigs.min5k_summary/{sample}.contigs.min5k_virus_summary.tsv",
        plasmid="output/MGE_detect/genomad_{sample}/{sample}.contigs.min5k_summary/{sample}.contigs.min5k_plasmid_summary.tsv"
    params:
        outdir="output/MGE_detect/genomad_{sample}",
        db="/lomi_home/zichun/MetaSpicii/db/genomad_db/"
    threads: config["threads"]
    conda:
        "MetaSpicii_version1"
    shell:
        """
        genomad end-to-end \
            {input} \
            {params.outdir} \
            {params.db} \
            --threads {threads}
        """
########################################
# 3️⃣ Extract virus + plasmid contigs
########################################

rule extract_mge_contigs:
    input:
        virus="output/MGE_detect/genomad_{sample}/{sample}.contigs.min5k_summary/{sample}.contigs.min5k_virus_summary.tsv",
        plasmid="output/MGE_detect/genomad_{sample}/{sample}.contigs.min5k_summary/{sample}.contigs.min5k_plasmid_summary.tsv"
    output:
        "output/MGE_detect/{sample}_MGE_contigs.tsv"
    shell:
        r"""
        (
            cut -f1 {input.virus} | tail -n +2
            cut -f1 {input.plasmid} | tail -n +2
        ) | awk '{{print $1"\t"NR}}' > {output}
        """

########################################
# 4️⃣ MetaTOR network construction
########################################

rule metator_network:
    input:
        assembly="output/MGE_detect/{sample}.contigs.min5k.fa",
        r1="qc/{sample}_hic_clean_1.fastq",
        r2="qc/{sample}_hic_clean_2.fastq"
    output:
        network="output/MGE_detect/metator_network_{sample}/network.txt",
        contigs="output/MGE_detect/metator_network_{sample}/contig_data_network.txt",
        pairs="output/MGE_detect/metator_network_{sample}/alignment_0_sorted.pairs.gz"
    params:
        outdir="output/MGE_detect/metator_network_{sample}"
    threads: config["threads"]
    conda:
        "MetaSpicii_version1"
    shell:
        """
        metator network \
          --assembly {input.assembly} \
          --forward {input.r1} \
          --reverse {input.r2} \
          --aligner bowtie2 \
          --aligner-mode normal \
          --outdir {params.outdir} \
          --threads {threads}
        """

########################################
# 5️⃣ MetaTOR partition
########################################

rule metator_partition:
    input:
        assembly="output/MGE_detect/{sample}.contigs.min5k.fa",
        contigs="output/MGE_detect/metator_network_{sample}/contig_data_network.txt",
        network="output/MGE_detect/metator_network_{sample}/network.txt"
    output:
        binning="output/MGE_detect/metator_partition_{sample}/binning.tsv"
    params:
        outdir="output/MGE_detect/metator_partition_{sample}"
    threads: config["threads"]
    conda:
        "MetaSpicii_version1"
    shell:
        """
        mkdir -p {params.outdir}

        # Step1: run metator partition
        metator partition \
          --assembly {input.assembly} \
          --contigs {input.contigs} \
          --network {input.network} \
          --outdir {params.outdir} \
          --threads {threads}

        # Step2: generate binning.tsv (FIXED)
        awk 'NR>1 && $8!="-" {{printf "%s\\t%d\\n", $2, int($8)}}' \
        {params.outdir}/contig_data_partition.txt > {output.binning}
        """
               
########################################
# 6️⃣ MetaTOR MGE linking
########################################

rule metator_mge_link:
    input:
        network="output/MGE_detect/metator_network_{sample}/network.txt",
        pairs="output/MGE_detect/metator_network_{sample}/alignment_0_sorted.pairs.gz",
        contigs_data="output/MGE_detect/metator_network_{sample}/contig_data_network.txt",
        binning="output/MGE_detect/metator_partition_{sample}/binning.tsv",
        mges="output/MGE_detect/{sample}_MGE_contigs.tsv",
        fasta="assembly/{sample}_final.contigs.fa"
    output:
        outdir=directory("output/MGE_detect/metator_mge_{sample}"),
        mge_fa="output/MGE_detect/metator_mge_{sample}/mges_binned.fa"
    threads: config["threads"]
    conda:
        "MetaSpicii_version1"
    shell:
        """
        metator mge \
          --network {input.network} \
          --binning {input.binning} \
          --mges {input.mges} \
          --contigs-data {input.contigs_data} \
          --fasta {input.fasta} \
          --outdir {output.outdir} \
          --threads {threads} \
          --threshold-bin 0.8 \
          --threshold-asso 0.1 \
          {input.pairs}
        """

