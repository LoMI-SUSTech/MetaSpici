# Assembly modes

MetaSpici assembly mode selects a route from the read inputs supplied to `metaspici run --mode assembly`. All routes use paired Hi-C reads and share the downstream reference association workflow.

| Route | Required reads | Main assembly/binning path | MGE discovery |
|---|---|---|---|
| Short reads + Hi-C | `--short-r1`, `--short-r2`, `--hic-r1`, `--hic-r2` | MEGAHIT; BASALT and MetaTOR binning; combine candidate bins; dRep dereplication | Contigs ≥5 kb → geNomad virus/plasmid identification → MetaTOR MGE linking |
| Long reads + Hi-C | `--long-reads`, `--long-type`, Hi-C pair | Flye or metaMDBG according to `--long-type`; MetaTOR binning; long-contig MAG candidates; dRep | Contigs ≥5 kb → geNomad → MetaTOR MGE linking |
| Hybrid + Hi-C | Short-read pair, long reads/type, Hi-C pair | OPERA-MS hybrid assembly with the associated long-read/merge workflow; combine MAG candidates; MetaTOR and dRep | Contigs ≥5 kb → geNomad → MetaTOR MGE linking |
| Hi-C only | Hi-C pair | HiCUP Truncater-derived read set, MEGAHIT, BASALT and MetaTOR binning; combine candidates; dRep | Contigs ≥5 kb → geNomad → MetaTOR MGE linking |

Supported `--long-type` values are `nanopore`, `pacbio`, `nanoporeR10`, and `pacbioHiFi`. Restriction enzyme labels accepted by the CLI are `Sau3AI`, `MboI`, and `DpnII`.

## Shared association stage

After assembly mode produces a MAG candidate set and MGE FASTA, it validates the FASTA inputs and invokes the shared reference workflow. The reference workflow performs Hi-C processing and mapping, ambiguity masking/filtering, contact reconstruction, depth filtering, fold-enrichment assessment, and final MGE–host table generation.

Reference-mode defaults are retained: minimum `Depth_Total` 0.15, strict fold enrichment `FE > 50`, and minimum contacts `>= 4`. Identity is specified as a percentage from 95 to 99.

## Important output paths

- `assembly_workflow/01_assembly/`: assembly products.
- `assembly_workflow/02_binning/`: individual and combined MAG candidates.
- `assembly_workflow/03_mge_detection/`: geNomad outputs and selected MGE candidates.
- `assembly_workflow/04_hic_network/`: MetaTOR network and partition outputs.
- `assembly_workflow/05_mge_host_linking/`: MGE binning/linking outputs.
- `assembly_workflow/02_binning/{sample}_final_binset/`: MAG candidates passed to the shared workflow.
- `assembly_workflow/05_mge_host_linking/metator_mge_{sample}/mges_binned.fa`: MGE FASTA passed to the shared workflow.
- `results/{sample}.associations.unfiltered_fe.tsv`: associations before FE filtering.
- `results/{sample}.mge_host.filtered.tsv`: final filtered MGE–host associations.

Use `metaspici run --help` for the current CLI arguments. Use a separate output directory for each sample and input configuration.
