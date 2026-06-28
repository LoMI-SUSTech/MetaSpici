configfile: "config.yaml"
wildcard_constraints:
    sample="[^/]+"
SAMPLES = [config["sample"]] 

rule all:
    input:
        expand("assembly/{sample}_final.contigs.fa", sample=SAMPLES),
        expand("binning/{sample}_final_shortbinset", sample=SAMPLES),
        expand("binning/{sample}_final_hicbinset", sample=SAMPLES),
        expand("binning/{sample}_final_binset", sample=SAMPLES),
        expand("output/MGE_detect/metator_mge_{sample}", sample=SAMPLES),
        expand("mapping/reference/reference_masked/{sample}_masked.fasta",sample=SAMPLES)

include: "../rules/hic_read_qc.py"
include: "../rules/hic_cutsite.py"
include: "../rules/short_assembly.py"
include: "../rules/length_filtering.py"
include: "../rules/short_binning.py"
include: "../rules/hic_binning.py"
include: "../rules/bin_drep.py"
include: "../rules/only_hic_mge.py"
include: "../rules/merged_links_short.py"
