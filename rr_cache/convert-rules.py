#!/usr/bin/env python3
import argparse
import gzip
import json
from collections import Counter
from typing import Dict, Any, List, Set, TextIO

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert metadata TSV(.gz) to (1) rr_reactions JSON and (2) template reactions JSON."
    )
    p.add_argument("-i", "--input", required=True, help="Input TSV file, optionally gzipped (.gz)")
    p.add_argument("--sep", default="\t", help="Field separator (default: tab)")
    p.add_argument("-o1", "--out-rr-reactions", required=True,
                   help="Output JSON for rr_reactions (TEMPLATE->RXN->entry); if ends with .gz, gzip it")
    p.add_argument("-o2", "--out-template-reactions", required=True,
                   help="Output JSON for template reactions (RXN->...); if ends with .gz, gzip it")
    p.add_argument("--only-valid", action="store_true", help="Keep only rows with VALID == True (if column exists)")
    return p.parse_args()


def open_text_auto(path: str, mode: str = "rt", encoding: str = "utf-8") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding=encoding)
    return open(path, mode, encoding=encoding)


def clean_str(x) -> str:
    if pd.isna(x):
        return ""
    return str(x)


def to_python_scalar(v):
    if pd.isna(v):
        return ""
    if hasattr(v, "item"):
        try:
            return v.item()
        except Exception:
            pass
    return v


def split_ids(dot_list: Any) -> List[str]:
    s = clean_str(dot_list).strip()
    if not s:
        return []
    return [p for p in s.split(".") if p]


def parse_side_counts(dot_list: Any) -> Dict[str, float]:
    parts = split_ids(dot_list)
    counts = Counter(parts)
    return {k: float(v) for k, v in counts.items()}


def direction_mnx(direction: str) -> int:
    d = (direction or "").strip().upper()
    if d in {"L2R", "LEFT_TO_RIGHT", "FORWARD"}:
        return 1
    if d in {"R2L", "RIGHT_TO_LEFT", "REVERSE"}:
        return -1
    return 0


def main():
    args = parse_args()

    # Read (plain or gz)
    with open_text_auto(args.input, "rt") as f:
        df = pd.read_csv(f, sep=args.sep)

    if args.only_valid and "VALID" in df.columns:
        df = df[df["VALID"] == True]

    # As requested:
    # - out-template-reactions: REACTION_ID -> {left,right,direction,main_left,main_right}
    # - out-rr-reactions: TEMPLATE_ID -> REACTION_ID -> detailed entry (+extra columns)
    template_reactions: Dict[str, Dict[str, Any]] = {}
    rr_reactions: Dict[str, Dict[str, Any]] = {}

    for _, row in df.iterrows():
        rxn_id = clean_str(row.get("REACTION_ID"))
        template_id = clean_str(row.get("TEMPLATE_ID"))

        left = parse_side_counts(row.get("LEFT_IDS"))
        right = parse_side_counts(row.get("RIGHT_IDS"))

        left_excl: Set[str] = set(split_ids(row.get("LEFT_EXCLUDED_IDS")))
        right_excl: Set[str] = set(split_ids(row.get("RIGHT_EXCLUDED_IDS")))

        main_left = [m for m in left.keys() if m not in left_excl]
        main_right = [m for m in right.keys() if m not in right_excl]

        # (2) template_reactions: per reaction
        template_reactions[rxn_id] = {
            "left": left,
            "right": right,
            "direction": direction_mnx(clean_str(row.get("DIRECTION"))),
            "main_left": main_left,
            "main_right": main_right,
        }

        # (1) rr_reactions: template -> rxn -> entry
        legacy = clean_str(row.get("LEGACY_ID"))
        subs_id = ""
        if "_" in legacy:
            subs_id = legacy.split("_")[-1].strip()
        if not subs_id:
            subs_id = next(iter(left.keys()), "")

        entry = {
            "rule_id": template_id,
            "rule_score": to_python_scalar(row.get("SCORE")),
            "reac_id": rxn_id,
            "subs_id": subs_id,
            "rel_direction": (1 if direction_mnx(clean_str(row.get("DIRECTION"))) >= 0 else -1),
            "left": left,
            "right": right,
        }

        # Add additional columns as new keys (lowercased), skipping those already mapped
        already_mapped = {"REACTION_ID", "TEMPLATE_ID", "SCORE", "LEFT_IDS", "RIGHT_IDS", "DIRECTION"}
        for col in df.columns:
            if col in already_mapped:
                continue
            key = col.lower()
            if key in entry:
                continue
            entry[key] = to_python_scalar(row.get(col))

        rr_reactions.setdefault(template_id, {})[rxn_id] = entry

    # Write single-line JSON (plain or gz)
    with open_text_auto(args.out_rr_reactions, "wt") as f:
        json.dump(rr_reactions, f, sort_keys=True, ensure_ascii=False)

    with open_text_auto(args.out_template_reactions, "wt") as f:
        json.dump(template_reactions, f, sort_keys=True, ensure_ascii=False)

    print(
        f"Wrote rr_reactions (TEMPLATE->RXN): {sum(len(v) for v in rr_reactions.values())} rows "
        f"across {len(rr_reactions)} templates -> {args.out_rr_reactions}\n"
        f"Wrote template_reactions (RXN->...): {len(template_reactions)} reactions -> {args.out_template_reactions}"
    )


if __name__ == "__main__":
    main()