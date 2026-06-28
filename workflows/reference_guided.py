configfile: "config.yaml"

wildcard_constraints:
    sample="[^./]+"

SAMPLES = [config["sample"]]


rule all:
    input:
        expand("hic_results/{sample}.final_result.tsv", sample=SAMPLES)


# ===== rules include（你后面自己补/替换）=====
include: "../rules/hic_reference_guided.py"
