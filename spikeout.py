#!/usr/bin/env python3
"""spikeout: quantify uniquely aligned reads across paired BAM files."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


@dataclass(frozen=True)
class BamDescriptor:
    path: Path
    sample_prefix: str
    assembly: str
    mapper: str


@dataclass(frozen=True)
class BamPair:
    sample_prefix: str
    mapper: str
    bam1: BamDescriptor
    bam2: BamDescriptor


def parse_bam_name(path: Path) -> BamDescriptor:
    stem = path.name
    if not stem.endswith(".bam"):
        raise ValueError(f"Input file is not a .bam file: {path}")

    name_without_ext = stem[:-4]
    parts = name_without_ext.rsplit("_", 2)
    if len(parts) != 3:
        raise ValueError(
            "Could not parse BAM name into '<sample_prefix>_<assembly>_<mapper>.bam': "
            f"{path.name}"
        )

    sample_prefix, assembly, mapper = parts
    if not sample_prefix or not assembly or not mapper:
        raise ValueError(f"Malformed BAM name: {path.name}")

    return BamDescriptor(path=path, sample_prefix=sample_prefix, assembly=assembly, mapper=mapper)


def build_pairs(paths: Iterable[Path]) -> list[BamPair]:
    grouped: dict[tuple[str, str], list[BamDescriptor]] = {}

    for path in paths:
        descriptor = parse_bam_name(path)
        grouped.setdefault((descriptor.sample_prefix, descriptor.mapper), []).append(descriptor)

    pairs: list[BamPair] = []
    errors: list[str] = []

    for (sample_prefix, mapper), descriptors in sorted(grouped.items()):
        if len(descriptors) != 2:
            errors.append(
                f"Expected exactly 2 BAMs for sample '{sample_prefix}' mapper '{mapper}', "
                f"found {len(descriptors)}"
            )
            continue

        sorted_descriptors = sorted(descriptors, key=lambda d: d.assembly)
        pairs.append(
            BamPair(
                sample_prefix=sample_prefix,
                mapper=mapper,
                bam1=sorted_descriptors[0],
                bam2=sorted_descriptors[1],
            )
        )

    if errors:
        raise ValueError("\n".join(errors))

    return pairs


def write_passing_read_names(bam_path: Path, mapq_cutoff: int, output_txt: Path) -> None:
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError(
            "pysam is required to parse BAM files. Install with: pip install pysam"
        ) from exc

    with pysam.AlignmentFile(str(bam_path), "rb") as bam, output_txt.open("w", encoding="utf-8") as out:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.mapping_quality < mapq_cutoff:
                continue
            out.write(f"{read.query_name}\n")


def sort_unique_input(input_txt: Path, output_txt: Path) -> None:
    sort_exe = shutil.which("sort")
    if not sort_exe:
        raise RuntimeError("GNU/BSD sort is required but was not found in PATH.")

    with output_txt.open("w", encoding="utf-8") as out_handle:
        subprocess.run(
            [sort_exe, "-u", str(input_txt)],
            check=True,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
        )


def iter_sorted_names(path: Path) -> Iterator[str]:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            yield line.rstrip("\n")


def count_set_differences(sorted_a: Path, sorted_b: Path) -> tuple[int, int]:
    a_iter = iter_sorted_names(sorted_a)
    b_iter = iter_sorted_names(sorted_b)

    unique_a = 0
    unique_b = 0

    a_name = next(a_iter, None)
    b_name = next(b_iter, None)

    while a_name is not None or b_name is not None:
        if b_name is None or (a_name is not None and a_name < b_name):
            unique_a += 1
            a_name = next(a_iter, None)
        elif a_name is None or b_name < a_name:
            unique_b += 1
            b_name = next(b_iter, None)
        else:
            a_name = next(a_iter, None)
            b_name = next(b_iter, None)

    return unique_a, unique_b


def quantify_pair(
    pair: BamPair,
    mapq_cutoff: int,
    temp_dir: Path,
    verbose: bool = False,
    pair_progress_label: str = "",
) -> tuple[int, int]:
    raw1 = temp_dir / f"{pair.sample_prefix}.{pair.bam1.assembly}.raw.txt"
    raw2 = temp_dir / f"{pair.sample_prefix}.{pair.bam2.assembly}.raw.txt"
    uniq1 = temp_dir / f"{pair.sample_prefix}.{pair.bam1.assembly}.uniq.txt"
    uniq2 = temp_dir / f"{pair.sample_prefix}.{pair.bam2.assembly}.uniq.txt"

    log_progress(
        verbose,
        f"{pair_progress_label} 1/4 extracting reads from {pair.bam1.path.name}",
    )
    write_passing_read_names(pair.bam1.path, mapq_cutoff, raw1)
    log_progress(
        verbose,
        f"{pair_progress_label} 2/4 extracting reads from {pair.bam2.path.name}",
    )
    write_passing_read_names(pair.bam2.path, mapq_cutoff, raw2)

    log_progress(verbose, f"{pair_progress_label} 3/4 sorting and deduplicating read names")
    sort_unique_input(raw1, uniq1)
    sort_unique_input(raw2, uniq2)

    raw1.unlink(missing_ok=True)
    raw2.unlink(missing_ok=True)

    log_progress(verbose, f"{pair_progress_label} 4/4 counting unique reads per assembly")
    return count_set_differences(uniq1, uniq2)


def validate_assemblies(pairs: list[BamPair]) -> tuple[str, str]:
    if not pairs:
        raise ValueError("No valid BAM pairs found.")

    first = (pairs[0].bam1.assembly, pairs[0].bam2.assembly)
    mismatched = [
        p for p in pairs if (p.bam1.assembly, p.bam2.assembly) != first
    ]
    if mismatched:
        raise ValueError(
            "Detected multiple assembly combinations across input pairs. "
            "This 3-column output format requires all pairs to share the same two assemblies."
        )

    return first


def write_output(
    pairs: list[BamPair],
    counts: list[tuple[int, int]],
    assemblies: tuple[str, str],
    output_handle: TextIO,
) -> None:
    output_handle.write(f"Sample\t{assemblies[0]}\t{assemblies[1]}\n")
    for pair, (count1, count2) in zip(pairs, counts):
        output_handle.write(f"{pair.sample_prefix}\t{count1}\t{count2}\n")


def log_progress(enabled: bool, message: str) -> None:
    if enabled:
        print(message, file=sys.stderr)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="spikeout",
        description=(
            "Quantify uniquely aligned reads between paired BAMs mapped to different assemblies."
        ),
    )
    parser.add_argument("bams", nargs="+", help="Input BAM files.")
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        default=20,
        help="MAPQ cutoff for considering a read aligned (default: 20).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output TSV path (default: stdout).",
    )
    parser.add_argument(
        "--tmpdir",
        default=".",
        help="Optional directory for temporary intermediate files (default current dir).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress and pairing summary messages.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    verbose = not args.quiet

    try:
        bam_paths = [Path(p).resolve() for p in args.bams]
        for p in bam_paths:
            if not p.exists():
                raise FileNotFoundError(f"Input BAM not found: {p}")

        pairs = build_pairs(bam_paths)
        assemblies = validate_assemblies(pairs)
        log_progress(
            verbose,
            (
                f"Found {len(pairs)} BAM pair(s) using assemblies: "
                f"{assemblies[0]} and {assemblies[1]}"
            ),
        )
        for i, pair in enumerate(pairs, start=1):
            log_progress(
                verbose,
                (
                    f"Pair {i}/{len(pairs)}: sample={pair.sample_prefix}, mapper={pair.mapper}, "
                    f"{pair.bam1.assembly}={pair.bam1.path.name}, "
                    f"{pair.bam2.assembly}={pair.bam2.path.name}"
                ),
            )

        with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
            temp_dir = Path(tmp)
            counts: list[tuple[int, int]] = []
            for i, pair in enumerate(pairs, start=1):
                log_progress(
                    verbose,
                    (
                        f"Processing pair {i}/{len(pairs)} ({pair.sample_prefix}): "
                        "extracting passing read names"
                    ),
                )
                count = quantify_pair(
                    pair,
                    args.mapq,
                    temp_dir,
                    verbose=verbose,
                    pair_progress_label=f"Pair {i}/{len(pairs)} ({pair.sample_prefix}):",
                )
                counts.append(count)
                log_progress(
                    verbose,
                    (
                        f"Processing pair {i}/{len(pairs)} ({pair.sample_prefix}): "
                        "complete"
                    ),
                )

        if args.output == "-":
            write_output(pairs, counts, assemblies, sys.stdout)
        else:
            out_path = Path(args.output)
            with out_path.open("w", encoding="utf-8") as out_handle:
                write_output(pairs, counts, assemblies, out_handle)

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
