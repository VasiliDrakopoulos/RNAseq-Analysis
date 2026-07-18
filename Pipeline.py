#!/usr/bin/env python3
"""
TRAP‑seq pipeline – full gene count matrix only.
NEBNext Low‑bias kit, paired‑end reads, mouse (GRCm38).
"""

import os, sys, subprocess, shutil
from pathlib import Path
import pandas as pd

def run_cmd(cmd, shell=False):
    print(f"  Running: {' '.join(cmd)}")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1, shell=shell)
    for line in proc.stdout:
        print("    " + line.rstrip())
    proc.wait()
    if proc.returncode != 0:
        print(f"  ERROR: command returned {proc.returncode}")
        sys.exit(1)

def check_tool(name):
    if shutil.which(name) is None:
        print(f"ERROR: {name} not found. Activate the conda environment first.")
        sys.exit(1)

# USER INPUTS
print("TRAP‑seq full gene pipeline (mouse)")
input_dir = input("Path to folder with sample subfolders: ").strip()
output_dir = input("Path where results will be saved: ").strip()

# Adapters (NEBNext Low‑bias kit)
r1_adapter_default = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
r2_adapter_default = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

a1 = input(f"Read 1 adapter [Enter for default: {r1_adapter_default[:30]}...]: ").strip()
if not a1:
    a1 = r1_adapter_default
a2 = input(f"Read 2 adapter [Enter for default: {r2_adapter_default[:30]}...]: ").strip()
if not a2:
    a2 = r2_adapter_default

hisat2_idx = input("HISAT2 index prefix (e.g., /home/vasili/mouse_hisat2_index/grcm38_tran/genome_tran): ").strip()
gtf = input("Annotation GTF (e.g., /home/vasili/mouse_reference/Mus_musculus.GRCm38.102.gtf): ").strip()

threads = input("Number of threads [4]: ").strip()
if not threads:
    threads = "4"

sample_input = input("Samples (name:group, space‑separated, e.g., sample1:A sample2:B): ").strip()

samples = []
for item in sample_input.split():
    if ':' in item:
        name, grp = item.split(':')
        samples.append((name.strip(), grp.strip()))
if not samples:
    print("No valid samples entered.")
    sys.exit(1)

input_path = Path(input_dir)
output_path = Path(output_dir)
output_path.mkdir(parents=True, exist_ok=True)

for tool in ['cutadapt', 'hisat2', 'samtools', 'featureCounts', 'fastqc', 'multiqc']:
    check_tool(tool)

# Find FASTQ files
sample_files = {}
for name, grp in samples:
    sample_dir = input_path / name
    r1 = list(sample_dir.glob("*_R1*.fastq.gz")) + list(sample_dir.glob("*_R1*.fq.gz"))
    r2 = list(sample_dir.glob("*_R2*.fastq.gz")) + list(sample_dir.glob("*_R2*.fq.gz"))
    if r1 and r2:
        sample_files[name] = (str(r1[0]), str(r2[0]))
    else:
        print(f"WARNING: {name} missing paired FASTQ files – skipping")

samples_final = [(name, grp) for name, grp in samples if name in sample_files]
if not samples_final:
    print("ERROR: No valid samples found.")
    sys.exit(1)

# OUTPUT DIRS
qc_dir = output_path / "qc"
trim_dir = output_path / "trimmed"
align_dir = output_path / "aligned"
counts_dir = output_path / "counts"
for d in [qc_dir, trim_dir, align_dir, counts_dir]:
    d.mkdir(exist_ok=True)

# FASTQC
print("\nStep 1: FastQC")
cmd = ['fastqc', '-o', str(qc_dir), '-t', threads]
for name, _ in samples_final:
    cmd.extend([sample_files[name][0], sample_files[name][1]])
run_cmd(cmd)
run_cmd(['multiqc', '-o', str(qc_dir), str(qc_dir)])

# CUTADAPT
print("\nStep 2: Adapter trimming")
trimmed = {}
for name, _ in samples_final:
    r1_in, r2_in = sample_files[name]
    r1_out = trim_dir / f"{name}_R1_trimmed.fastq.gz"
    r2_out = trim_dir / f"{name}_R2_trimmed.fastq.gz"
    run_cmd(['cutadapt', '-a', a1, '-A', a2, '-o', str(r1_out), '-p', str(r2_out),
             '-m', '15', '-q', '20', r1_in, r2_in])
    trimmed[name] = (str(r1_out), str(r2_out))

# HISAT2
print("\nStep 3: HISAT2 alignment")
for name, _ in samples_final:
    r1, r2 = trimmed[name]
    sorted_bam = align_dir / f"{name}_sorted.bam"
    log_file = align_dir / f"{name}_hisat2.log"
    hisat_cmd = ['hisat2', '-x', hisat2_idx, '-1', r1, '-2', r2,
                 '-p', threads, '--dta', '--new-summary']
    print(f"  Aligning {name}")
    with open(log_file, 'w') as log_f:
        p1 = subprocess.Popen(hisat_cmd, stdout=subprocess.PIPE, stderr=log_f)
        p2 = subprocess.Popen(['samtools', 'sort', '-@', threads, '-o', str(sorted_bam), '-'],
                              stdin=p1.stdout, stderr=subprocess.PIPE)
        p1.stdout.close()
        p2.communicate()
        p1.wait()
    run_cmd(['samtools', 'index', str(sorted_bam)])

# FEATURECOUNTS (inclusive)
print("\nStep 4: Gene counting (featureCounts)")
bams = [str(align_dir / f"{name}_sorted.bam") for name, _ in samples_final]
counts_file = counts_dir / "raw_counts.txt"
run_cmd(['featureCounts', '-p', '-t', 'exon', '-g', 'gene_id',
         '-M', '--fraction', '-O', '-a', gtf, '-o', str(counts_file), *bams])

# BUILD FULL GENE COUNT MATRIX
print("\nStep 5: Creating full gene count matrix")
df = pd.read_csv(counts_file, sep='\t', comment='#', header=0)
count_cols = df.columns[6:]
mat = df.set_index('Geneid')[count_cols]
mat.columns = [c.split('/')[-1].replace('_sorted.bam', '') for c in mat.columns]
mat = mat[mat.sum(axis=1) > 0]
mat.to_csv(output_path / "full_gene_counts.csv")
print(f"  Full gene counts saved: {output_path / 'full_gene_counts.csv'}")
print(f"  Number of genes with >0 reads: {mat.shape[0]}")

print("\n===== PIPELINE FINISHED =====")
print(f"Results are in: {output_path}")