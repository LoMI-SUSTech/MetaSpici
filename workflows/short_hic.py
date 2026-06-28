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
        expand("merged_links.tsv",sample=SAMPLES ),
        expand("hic_results/{sample}.final_result.tsv",sample=SAMPLES)


include: "../rules/short_read_qc.py"
include: "../rules/short_assembly.py"
include: "../rules/length_filtering.py"
include: "../rules/short_binning.py"
include: "../rules/hic_read_qc.py"  
include: "../rules/hic_binning.py"
include: "../rules/bin_drep_version1.py"
include: "../rules/only_hic_mge.py"
include: "../rules/merged_links_short_version2.py"

