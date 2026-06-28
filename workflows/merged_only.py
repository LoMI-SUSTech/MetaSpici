configfile: "config.yaml"

wildcard_constraints:
    sample="[^/]+"

SAMPLES = [config["sample"]]

# 只加载你需要的规则
include: "../rules/merged_links_short_version1.py"

# 只定义最终目标（不会牵连其他模块）
rule all:
    input:
        expand(
            "merged_links.tsv",
            sample=SAMPLES
        )
