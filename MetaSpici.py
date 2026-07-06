#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
import yaml


def get_script_dir():
    return os.path.dirname(os.path.abspath(sys.argv[0]))


def detect_mode(args):
    # ⭐ reference-guided 模式（原 user-mode）
    if args.reference_guided:
        if not (args.mag_dir and args.votu):
            print("Error: --mag_dir and --votu are required in --reference-guided mode")
            sys.exit(1)
        if not (args.hic_read1 and args.hic_read2):
            print("Error: HiC data is required in reference-guided mode")
            sys.exit(1)
        return "reference_guided", None, None

    # ===== 原有逻辑（完全保留）=====
    if not (args.hic_read1 and args.hic_read2):
        print("Error: HiC data is required. Please provide both --hic-read1 and --hic-read2")
        sys.exit(1)

    has_short = bool(args.read1 and args.read2)
    long_read = None
    long_type = None

    if args.nanopore:
        long_read = args.nanopore
        long_type = "nanopore"
    elif args.nanopore_R10:
        long_read = args.nanopore_R10
        long_type = "nanoporeR10"
    elif args.pacbio:
        long_read = args.pacbio
        long_type = "pacbio"
    elif args.pacbioHiFi:
        long_read = args.pacbioHiFi
        long_type = "pacbioHiFi"

    has_long = bool(long_read)

    if has_short and has_long:
        mode_str = "short_long_hic"
    elif has_short and not has_long:
        mode_str = "short_hic"
    elif not has_short and has_long:
        mode_str = "long_hic"
    elif not has_short and not has_long:
        mode_str = "hic_only_mge"
    else:
        mode_str = "hic_only_mge"

    return mode_str, long_read, long_type


def main():
    parser = argparse.ArgumentParser(
        description="MetaSpicii - Multi-modal Genome Assembly and Hi-C Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # ===== short reads =====
    parser.add_argument("-s1", "--read1", help="Short Read1 file")
    parser.add_argument("-s2", "--read2", help="Short Read2 file")

    # ===== long reads =====
    long_group = parser.add_mutually_exclusive_group()
    long_group.add_argument("--nanopore", help="ONT Nanopore reads")
    long_group.add_argument("--nanopore_R10", help="ONT Nanopore reads R10.4+")
    long_group.add_argument("--pacbioHiFi", help="PacBio HiFi reads")
    long_group.add_argument("--pacbio", help="PacBio CLR reads")

    # ===== Hi-C =====
    parser.add_argument("-h1", "--hic_read1", required=True, help="Hi-C read R1 (required)")
    parser.add_argument("-h2", "--hic_read2", required=True, help="Hi-C read R2 (required)")
    parser.add_argument("-e", "--enzyme", required=True, help="Hi-C digest enzyme")

    # ===== ⭐ reference-guided 模式 =====
    parser.add_argument("--reference-guided", action="store_true",
                        help="Use reference-guided mode (MAGs + vOTUs + Hi-C)")
    parser.add_argument("--mag_dir", help="Directory containing MAGs")
    parser.add_argument("--votu", help="vOTU fasta file")

    # ===== 基础参数 =====
    parser.add_argument("-n", "--sample", required=True, help="Sample name")
    parser.add_argument("-j", "--threads", default="1", help="Number of threads")
    parser.add_argument("-o", "--output", default=".", help="Output directory")
    parser.add_argument("--target", help="Snakemake target rule (optional)")
    parser.add_argument("--create-envs-only", action="store_true",
                        help="Only create conda environments")

    args = parser.parse_args()

    mode, long_read, long_type = detect_mode(args)

    os.makedirs(args.output, exist_ok=True)

    # ===== config.yaml =====
    config_data = {
        "sample": args.sample,
        "mode": mode,
        "threads": int(args.threads),
        "output_dir": args.output,
        "enzyme": args.enzyme,
    }

    # ===== reference-guided 输入 =====
    if args.reference_guided:
        config_data["reference_inputs"] = {
            "mag_dir": os.path.abspath(args.mag_dir),
            "votu": os.path.abspath(args.votu)
        }

    # ===== short reads =====
    if args.read1 and args.read2:
        config_data["short_reads"] = {
            "r1": os.path.abspath(args.read1),
            "r2": os.path.abspath(args.read2)
        }

    # ===== long reads =====
    if long_read:
        config_data["long_reads"] = {
            "type": long_type,
            "file": os.path.abspath(long_read)
        }

    # ===== Hi-C =====
    config_data["hic_reads"] = {
        "r1": os.path.abspath(args.hic_read1),
        "r2": os.path.abspath(args.hic_read2)
    }

    # ===== workflow =====
    script_dir = get_script_dir()
    workflows_dir = os.path.join(script_dir, "workflows")
    workflow_file = os.path.join(workflows_dir, f"{mode}.py")

    config_file = os.path.join(args.output, "config.yaml")

    with open(config_file, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

    # ===== snakemake =====
    cmd = [
        "snakemake",
        "-s", workflow_file,
        "--cores", str(args.threads),
        "--directory", os.path.abspath(args.output),
        "--configfile", os.path.abspath(config_file),
        "--use-conda",
    ]

    if args.target:
        cmd.append(args.target)

    if args.create_envs_only:
        cmd.append("--conda-create-envs-only")

    print(f"\n{'='*60}")
    print(f"Running MetaSpicii in {mode} mode")
    print(f"Sample: {args.sample}")
    print(f"Threads: {args.threads}")
    print(f"Output directory: {args.output}")
    print(f"Long read type: {long_type if long_type else 'None'}")

    if args.reference_guided:
        print(f"MAG directory: {args.mag_dir}")
        print(f"vOTU file: {args.votu}")

    print(f"{'='*60}\n")
    print(f"Command: {' '.join(cmd)}\n")

    subprocess.run(cmd)


if __name__ == "__main__":
    main()
