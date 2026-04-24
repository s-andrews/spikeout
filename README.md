# spikeout

A memory-efficient command-line program to quantify uniquely aligned reads across paired BAM files from spike-in experiments.

## Input naming convention

`spikeout` expects BAM names in this format:

`<sample_prefix>_<assembly>_<mapper>.bam`

Example pair:

- `EZH2_Y726A_EZH1_KO_H3K27me3_rep3_dec19_BDGP6_bowtie2.bam`
- `EZH2_Y726A_EZH1_KO_H3K27me3_rep3_dec19_GRCm39_bowtie2.bam`

The program matches files using the same `<sample_prefix>` and `<mapper>`, then compares the two assemblies.

## Usage

```bash
python spikeout.py *.bam --mapq 20 --output spikeout_counts.tsv
```

Options:

- `--mapq / -q`: MAPQ cutoff (default: `20`)
- `--output / -o`: output TSV path (default: stdout)
- `--tmpdir`: optional temp directory for intermediates

## Output

A tab-delimited file with 3 columns:

1. `Sample`
2. unique aligned read count for assembly 1
3. unique aligned read count for assembly 2

Header assemblies are inferred from BAM filenames (alphabetical order).

## Requirements

- Python 3.9+
- `pysam`
- `sort` available on `PATH`
