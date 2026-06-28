configfile: "config.yaml"
SAMPLES = [config["sample"]] 

rule all:
    input:
        expand("assembly/{sample}_final.contigs.fa", sample=SAMPLES),
        expand("binning/{sample}_final_shortbinset", sample=SAMPLES),
        expand("binning/{sample}_final_hicbinset", sample=SAMPLES),
        expand("binning/{sample}_final_binset", sample=SAMPLES)

include: "../rules/hic_read_qc.py"
include: "../rules/hic_cutsite.py"
include: "../rules/short_assembly.py"
include: "../rules/length_filtering.py"
include: "../rules/short_binning.py"
include: "../rules/hic_binning.py"
include: "../rules/bin_drep.py"