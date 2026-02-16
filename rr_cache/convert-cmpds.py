#!/usr/bin/env python3
import argparse
import gzip
import io
import json
from typing import TextIO

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Convert compounds TSV(.gz) to JSON for VALID chemicals only.")
    p.add_argument("-i", "--input", required=True, help="Input TSV file, optionally gzipped (.gz)")
    p.add_argument("-o", "--output", required=True, help="Output JSON file; if ends with .gz, gzip it")
    p.add_argument("--sep", default="\t", help="Field separator (default: tab)")
    p.add_argument("--only-valid", action="store_true", help="Keep only rows with VALID == True")
    return p.parse_args()


def open_text_auto(path: str, mode: str = "rt", encoding: str = "utf-8") -> TextIO:
    """
    Open plain or gzipped text based on filename suffix.
    mode: 'rt' for read text, 'wt' for write text.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding=encoding)
    return open(path, mode, encoding=encoding)


def main():
    args = parse_args()

    # Read (plain or gz)
    with open_text_auto(args.input, "rt") as f:
        df = pd.read_csv(f, sep=args.sep)

    if args.only_valid and "VALID" in df.columns:
        df = df[df["VALID"] == True]

    output = {}
    for _, row in df.iterrows():
        mnx_id = row["ID"]

        smiles = "" if pd.isna(row.get("SMILES")) else str(row.get("SMILES"))
        smiles_0 = "" if pd.isna(row.get("SMILES_0")) else str(row.get("SMILES_0"))

        entry = {"smiles": smiles}
        if smiles_0 and smiles_0 != smiles:
            entry["smiles_0"] = smiles_0

        output[str(mnx_id)] = entry

    # Write (plain or gz) single-line JSON
    with open_text_auto(args.output, "wt") as f:
        json.dump(output, f, sort_keys=True, ensure_ascii=False)

    print(f"Wrote {len(output)} chemicals to {args.output}")


if __name__ == "__main__":
    main()