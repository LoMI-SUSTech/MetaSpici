configfile: "config.yaml"
SAMPLES = [config["sample"]] 

rule all:
    input:
        expand("assembly/{sample}_final.contigs.fa", sample=SAMPLES),
        expand("assembly/{sample}_single_genome/split_contigs", sample=SAMPLES),
        expand("binning/{sample}_final_hicbinset", sample=SAMPLES),
        expand("binning/{sample}_final_binset", sample=SAMPLES)

include: "../rules/long_assembly.py"
include: "../rules/length_filtering.py"
include: "../rules/getSingle_genome.py"
include: "../rules/hic_binning.py"
include: "../rules/hic_read_qc.py"
include: "../rules/bin_drep.py"