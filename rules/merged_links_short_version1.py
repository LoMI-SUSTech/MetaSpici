import glob

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

rule fastuniq:
    input:
        r1="qc/{sample}_hic_clean_1.fastq",
        r2="qc/{sample}_hic_clean_2.fastq"
    output:
        r1="mapping/fastuniq/{sample}_R1.fastq",
        r2="mapping/fastuniq/{sample}_R2.fastq"
    conda:
        "fastuniq"
    shell:
        """
        mkdir -p mapping/fastuniq

        echo {input.r1} > mapping/fastuniq/{wildcards.sample}_list.txt
        echo {input.r2} >> mapping/fastuniq/{wildcards.sample}_list.txt

        fastuniq \
            -i mapping/fastuniq/{wildcards.sample}_list.txt \
            -t q \
            -o {output.r1} \
            -p {output.r2}

        rm mapping/fastuniq/{wildcards.sample}_list.txt
        """

rule hic_cutsite_trim:
    input:
        r1="mapping/fastuniq/{sample}_R1.fastq",
        r2="mapping/fastuniq/{sample}_R2.fastq"
    output:
        r1="mapping/cutsite/{sample}_R1.fq.gz",
        r2="mapping/cutsite/{sample}_R2.fq.gz"
    conda:
        "hicstuff"
    shell:
        """
        mkdir -p mapping/cutsite
        hicstuff cutsite -1 {input.r1} -2 {input.r2} -e Sau3AI -t 16 -p mapping/cutsite/{wildcards.sample}
        """

rule bwa_index:
    input:
        "mapping/reference/reference_masked/{sample}_masked.fasta"
    output:
        "mapping/reference/reference_masked/{sample}_masked.fasta.bwt"
    conda:
        "map"
    shell:
        """
        bwa index {input}
        """

rule map_reads:
    input:
        r1="mapping/cutsite/{sample}_R1.fq.gz",
        r2="mapping/cutsite/{sample}_R2.fq.gz",
        ref="mapping/reference/reference_masked/{sample}_masked.fasta",
        idx="mapping/reference/reference_masked/{sample}_masked.fasta.bwt"
    output:
        bed="mapped_reads/{sample}.bed"
    threads: config["threads"]
    conda:
        "map"
    shell:
        """
        mkdir -p mapped_reads
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samblaster | \
        samtools view -S -h -b -q 20 -F 2316 | \
        samtools sort -o mapped_reads/{wildcards.sample}.sorted.bam
        bedtools bamtobed -i mapped_reads/{wildcards.sample}.sorted.bam > {output.bed}
        """

rule split_fasta:
    input:
        "mapping/reference/reference_masked/{sample}_masked.fasta"
    output:
        directory("blast_input_mock/{sample}")
    conda:"seqkit"
    shell:
        """
        mkdir -p {output}
        seqkit split2 -s 1 -O {output} {input}
        """

#######################################
# 7️⃣ Pairwise BLAST of single-sequence FASTAs
#######################################
rule pairwise_blast:
    input:
        expand("blast_input_mock/{sample}", sample=SAMPLES)
    output:
        directory("blast_results_mock")
    threads: config["threads"]
    conda: "blast"
    shell:
        """
        mkdir -p blast_results_mock
        FILES=(blast_input_mock/*/*.fasta)

        for ((i=0; i<${{#FILES[@]}}; i++)); do
          for ((j=i+1; j<${{#FILES[@]}}; j++)); do

            q=${{FILES[$i]}}
            db=${{FILES[$j]}}

            qname=$(basename $q .fasta)
            dbname=$(basename $db .fasta)

            makeblastdb -in $db -dbtype nucl -out $db.db

            blastn -query $q -db $db.db \
              -out blast_results_mock/${{qname}}_vs_${{dbname}}.tsv \
              -outfmt '6 qseqid sseqid pident length qstart qend sstart send' \
              -evalue 1e-10 \
              -num_threads {threads}

          done
        done
        """

#######################################
# 8️⃣ Merge BLAST regions
#######################################
rule merge_blast_regions:
    input:
        "blast_results_mock"
    output:
        "merged_regions_mock.bed"
    conda: "map"
    shell:
        r"""
        cat blast_results_mock/*.tsv > merged_filtered_blast_mock.tsv.tmp

        awk '($1 != $2) && ($3 > 99) && ($4 > 75)' \
            merged_filtered_blast_mock.tsv.tmp > merged_filtered_blast_mock.tsv

        awk -F'\t' '{{s=$5;e=$6;if(s>e){{tmp=s;s=e;e=tmp}}; print $1"\t"s"\t"e;
                     s=$7;e=$8;if(s>e){{tmp=s;s=e;e=tmp}}; print $2"\t"s"\t"e}}' \
            merged_filtered_blast_mock.tsv > 3column_mock.tsv

        sort -k1,1 -k2,2n 3column_mock.tsv > 3column_mock.sorted.tsv

        bedtools merge -i 3column_mock.sorted.tsv > {output}

        rm merged_filtered_blast_mock.tsv.tmp 3column_mock.tsv 3column_mock.sorted.tsv
        """

#######################################
# 9️⃣ Filter mapping BEDs with repeats
#######################################
rule filter_mapping_bed:
    input:
        bed=lambda wildcards: f"mapped_reads/{wildcards.sample}.bed",
        repeat="merged_regions_mock.bed"
    output:
        "filtered_beds_clean/{sample}.filtered.bed"
    conda: "map"
    shell:
        """
        mkdir -p filtered_beds_clean
        bedtools subtract \
            -a {input.bed} \
            -b {input.repeat} \
            -nonamecheck \
            > {output}
        """

rule get_circo:
    input:
        "filtered_beds_clean/{sample}.filtered.bed"
    output:
        diff="filtered_beds_clean/{sample}.filtered.bed.ciro_diff.txt",
        same="filtered_beds_clean/{sample}.filtered.bed.ciro_same.txt"
    conda:"openpyx"
    shell:
        """
        python scripts/get_circo_file.py {input}
        """

rule deduplicate_diff:
    input:
        "filtered_beds_clean/{sample}.filtered.bed.ciro_diff.txt"
    output:
        "filtered_beds_clean/{sample}.filtered.bed.ciro_diff.txt.reads_x_cleaned.txt"
    conda:"openpyx"
    shell:
        """
        python scripts/deduplicate.py {input}
        """

rule deduplicate_same:
    input:
        "filtered_beds_clean/{sample}.filtered.bed.ciro_same.txt"
    output:
        "filtered_beds_clean/{sample}.filtered.bed.ciro_same.txt.reads_x_cleaned.txt"
    conda:"openpyx"
    shell:
        """
        python scripts/deduplicate.py {input}
        """

rule merge_ciro_files:
    input:
        diff="filtered_beds_clean/{sample}.filtered.bed.ciro_diff.txt.reads_x_cleaned.txt",
        same="filtered_beds_clean/{sample}.filtered.bed.ciro_same.txt.reads_x_cleaned.txt"
    output:
        "filtered_beds_clean/{sample}.bed_total_duplicates_modified"
    shell:
        """
        cat {input.diff} {input.same} \
        > filtered_beds_clean/{wildcards.sample}.bed_total_duplicates

        sed 's/$/:00/' \
        filtered_beds_clean/{wildcards.sample}.bed_total_duplicates \
        > {output}
        """

rule build_contig_table:
    input:
        bed="filtered_beds_clean/{sample}.filtered.bed"
    output:
        "filtered_beds_clean/{sample}.merged_output.txt"
    shell:
        """
        awk -v FS='\\t' -v OFS='\\t' '{print $1,$4}' {input.bed} | \
        awk -F'/' '{print $1}' | \
        awk -v FS='\\t' -v OFS='\\t' '{print $2,$1}' > filtered_beds_clean/{wildcards.sample}.aligned.contig

        awk '{count[$1]++} END {for (id in count) print id, count[id]}' \
        filtered_beds_clean/{wildcards.sample}.aligned.contig \
        > filtered_beds_clean/{wildcards.sample}.read_counts.txt

        awk '$2 > 1' \
        filtered_beds_clean/{wildcards.sample}.read_counts.txt \
        > filtered_beds_clean/{wildcards.sample}.filtered_reads.txt

        awk 'NR==FNR{a[$1];next} $1 in a' \
        filtered_beds_clean/{wildcards.sample}.filtered_reads.txt \
        filtered_beds_clean/{wildcards.sample}.aligned.contig \
        > {output}
        """

rule build_final_output:
    input:
        "filtered_beds_clean/{sample}.merged_output.txt"
    output:
        "filtered_beds_clean/{sample}.final_output.txt"
    shell:
        """
        awk -v OFS='\\t' '
        {
            read_id = $1
            contig_id = $2
            contigs[read_id][count[read_id]++] = contig_id
        }
        END {
            for (id in contigs) {
                if (count[id] >= 2) {
                    asort(contigs[id])
                    print id, contigs[id][1], contigs[id][2], "1"
                }
            }
        }
        ' {input} > {output}
        """

rule filter_final_output:
    input:
        final="filtered_beds_clean/{sample}.final_output.txt",
        dup="filtered_beds_clean/{sample}.bed_total_duplicates_modified"
    output:
        "filtered_beds_clean/{sample}.final_output_filtered.txt"
    shell:
        """
        awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' \
        {input.dup} {input.final} > {output}
        """

rule merge_links:
    input:
        expand("filtered_beds_clean/{sample}.final_output_filtered.txt", sample=SAMPLES)
    output:
        "merged_links.tsv"
    shell:
        """
        awk -F '\\t' '{key=$2"\\t"$3; sum[key]+=$4} END {for (k in sum) print k"\\t"sum[k]}' \
        filtered_beds_clean/*.final_output_filtered.txt > {output}
        """
