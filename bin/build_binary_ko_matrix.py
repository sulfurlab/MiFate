#!/usr/bin/env python3
import sys
import csv
import gzip
import re
import argparse
import os

KO_PAT = re.compile(r"^K\d+$")


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def normalize_mag_name(x):
    """
    把输入统一标准化为 MAG 名：
    /path/to/A.fna -> A
    A.fna -> A
    A -> A
    """
    x = x.strip()
    x = os.path.basename(x)
    x = re.sub(r"\.(fna|fa|fasta|faa|gz)$", "", x)
    return x


def read_mag_list(mag_list_path):
    target_order = []
    target_set = set()

    with open(mag_list_path, "r") as f:
        for line in f:
            raw = line.strip()
            if not raw:
                continue

            g = normalize_mag_name(raw)
            if not g:
                continue

            if g not in target_set:
                target_set.add(g)
                target_order.append(g)

    return target_order, target_set


def parse_mapper(mapper_path, target_set):
    rows = {g: set() for g in target_set}
    all_kos = set()
    matched_genomes = set()

    with open_any(mapper_path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue

            parts = ln.split("\t")
            left = parts[0]
            ko_field = parts[1] if len(parts) > 1 else ""

            genome = left.split("|", 1)[0].strip()
            if genome not in target_set:
                continue

            matched_genomes.add(genome)

            if ko_field:
                for k in re.split(r"[ ,;]+", ko_field.strip()):
                    if KO_PAT.match(k):
                        rows[genome].add(k)
                        all_kos.add(k)

    return rows, sorted(all_kos), matched_genomes


def write_binary_matrix(out_csv, target_order, rows, kos):
    with open(out_csv, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["Genome"] + kos)
        for g in target_order:
            presence = [1 if k in rows.get(g, set()) else 0 for k in kos]
            w.writerow([g] + presence)


def main():
    parser = argparse.ArgumentParser(
        description="Build Genome x KO binary matrix (0/1) from kofam mapper."
    )
    parser.add_argument("-m", "--mapper", required=True, help="Input kofam.mapper or .gz")
    parser.add_argument("-l", "--mag_list", required=True, help="MAG list file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    target_order, target_set = read_mag_list(args.mag_list)
    if not target_order:
        sys.stderr.write("[ERR] MAG list is empty after normalization.\n")
        sys.exit(1)

    rows, kos, matched_genomes = parse_mapper(args.mapper, target_set)

    if not matched_genomes:
        sys.stderr.write("[ERR] No MAG names in mag_list matched any genome in mapper.\n")
        sys.stderr.write("[TIP] Check whether mag_list contains paths or filenames inconsistent with mapper headers.\n")
        sys.exit(2)

    if not kos:
        sys.stderr.write("[WARN] MAGs were matched, but no valid KO IDs were found.\n")

    write_binary_matrix(args.output, target_order, rows, kos)

    print(f"[DONE] Binary matrix written to: {args.output}")
    print(f"[INFO] Number of MAGs in list: {len(target_order)}")
    print(f"[INFO] Number of matched MAGs: {len(matched_genomes)}")
    print(f"[INFO] Number of KOs: {len(kos)}")


if __name__ == "__main__":
    main()