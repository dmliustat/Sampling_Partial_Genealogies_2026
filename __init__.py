from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .rank_distance import (
    RankDistanceAnalysis,
    analyze_tree_sequence_folder,
    discover_tree_sequence_files,
    rank_distance,
    simulate_reference_tree_sequence,
    summarize_tree_sequence,
)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="dataanalysis",
        description="Analyze partial .trees files with rank-based summaries.",
    )
    parser.add_argument(
        "folder",
        nargs="?",
        default="InputFiles",
        help="Folder containing .trees files or partial tree-sequence outputs.",
    )
    parser.add_argument("--n", type=int, required=True, help="Global lineage count for every file.")
    parser.add_argument(
        "--n-tau",
        type=int,
        required=True,
        dest="n_tau",
        help="Global lineage count at the stopping time tau for every file.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search subdirectories when discovering input files.",
    )
    args = parser.parse_args(argv)

    result = analyze_tree_sequence_folder(
        Path(args.folder),
        recursive=args.recursive,
        sample_count_config=args.n,
        stopping_lineages_config=args.n_tau,
    )
    print(result.format_matrix())
    return 0


__all__ = [
    "RankDistanceAnalysis",
    "analyze_tree_sequence_folder",
    "discover_tree_sequence_files",
    "main",
    "rank_distance",
    "simulate_reference_tree_sequence",
    "summarize_tree_sequence",
]
