from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import msprime
import tskit


@dataclass(frozen=True)
class RankSegmentSummary:
    left: float
    right: float
    tree_index: int
    span: float
    root_count: int
    observed_rank_count: int
    missing_rank_count: int
    observed_rank_mass: int
    missing_rank_mass: int
    stopping_time: float | None


@dataclass(frozen=True)
class TreeSequenceRankSummary:
    path: Path | None
    sample_count_config: int
    stopping_lineages_config: int
    sample_count: int
    sequence_length: float
    tree_count: int
    full_rank_mass: int
    observed_rank_mass_expected: int
    ancient_rank_mass_expected: int
    ancient_rank_value: float
    segments: tuple[RankSegmentSummary, ...]
    mean_root_count: float
    mean_observed_rank_count: float
    mean_stopping_time: float | None
    observed_rank_mass_fraction: float
    missing_rank_mass_fraction: float


@dataclass(frozen=True)
class RankDistanceAnalysis:
    paths: tuple[Path, ...]
    summaries: tuple[TreeSequenceRankSummary, ...]
    matrix: tuple[tuple[float, ...], ...]

    def format_matrix(self) -> str:
        labels = [path.name for path in self.paths]
        if not labels:
            return "No tree-sequence files were found."

        header = ["file"] + labels
        rows = ["\t".join(header)]
        for label, row in zip(labels, self.matrix, strict=True):
            formatted = [f"{value:.6f}" for value in row]
            rows.append("\t".join([label] + formatted))
        return "\n".join(rows)


def discover_tree_sequence_files(
    folder: Path | str,
    recursive: bool = False,
    suffixes: tuple[str, ...] = (".trees",),
) -> tuple[Path, ...]:
    root = Path(folder)
    if not root.exists():
        raise FileNotFoundError(f"Input folder does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Input path is not a folder: {root}")

    candidates = root.rglob("*") if recursive else root.iterdir()
    files = sorted(
        path
        for path in candidates
        if path.is_file()
        and not path.name.startswith(".")
        and (not suffixes or path.suffixes and any(path.name.endswith(suffix) for suffix in suffixes))
    )
    return tuple(files)


def load_tree_sequence(path: Path | str) -> tskit.TreeSequence:
    source = Path(path)
    try:
        return tskit.load(source)
    except Exception as error:  # pragma: no cover - re-raised with context
        raise ValueError(f"Could not load tree sequence from {source}") from error


def simulate_reference_tree_sequence(
    sample_size: int = 4,
    sequence_length: float = 1.0,
    recombination_rate: float = 0.0,
    population_size: float = 10_000,
    random_seed: int = 12345,
) -> tskit.TreeSequence:
    return msprime.sim_ancestry(
        samples=sample_size,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        population_size=population_size,
        random_seed=random_seed,
    )


def ancient_rank_value(sample_count: int, stopping_lineages: int) -> float:
    if sample_count < 2:
        raise ValueError("sample_count must be at least 2")
    if stopping_lineages < 2:
        raise ValueError("stopping_lineages must be at least 2")
    if stopping_lineages > sample_count:
        raise ValueError("stopping_lineages cannot exceed sample_count")

    total_rank_mass = sample_count * (sample_count - 1) / 2
    observed_rank_mass = (sample_count - stopping_lineages) * (sample_count - stopping_lineages + 1) / 2
    ancient_count = stopping_lineages - 1
    return (total_rank_mass - observed_rank_mass) / ancient_count


def summarize_tree_sequence(
    path_or_ts: Path | str | tskit.TreeSequence,
    /,
    sample_count_config: int | None = None,
    stopping_lineages_config: int | None = None,
) -> TreeSequenceRankSummary:
    if isinstance(path_or_ts, tskit.TreeSequence):
        tree_sequence = path_or_ts
        source_path: Path | None = None
    else:
        source_path = Path(path_or_ts)
        tree_sequence = load_tree_sequence(source_path)

    sample_count = tree_sequence.num_samples
    if sample_count_config is None:
        sample_count_config = sample_count
    if stopping_lineages_config is None:
        if tree_count == 0:
            raise ValueError("Tree sequence must contain at least one tree")
        stopping_lineages_config = max(len(next(tree_sequence.trees()).roots), 2)
    if sample_count != sample_count_config:
        raise ValueError(
            f"Tree sequence sample count {sample_count} does not match configured n={sample_count_config}"
        )
    if stopping_lineages_config < 2:
        raise ValueError("stopping_lineages_config must be at least 2")
    if stopping_lineages_config > sample_count_config:
        raise ValueError("stopping_lineages_config cannot exceed sample_count_config")

    tree_count = tree_sequence.num_trees
    full_rank_mass = sample_count * (sample_count - 1) // 2
    observed_rank_mass_expected = (sample_count - stopping_lineages_config) * (sample_count - stopping_lineages_config + 1) // 2
    ancient_rank_mass_expected = full_rank_mass - observed_rank_mass_expected
    ancient_rank_value_expected = ancient_rank_value(sample_count, stopping_lineages_config)

    segments: list[RankSegmentSummary] = []
    weighted_root_count = 0.0
    weighted_observed_rank_count = 0.0
    weighted_observed_rank_mass = 0.0
    weighted_missing_rank_mass = 0.0
    weighted_stopping_time = 0.0
    stopping_time_is_defined = False

    for tree_index, tree in enumerate(tree_sequence.trees()):
        root_count = len(tree.roots)
        observed_rank_count = max(sample_count - root_count, 0)
        missing_rank_count = max(root_count - 1, 0)
        observed_rank_mass = observed_rank_mass_expected
        missing_rank_mass = ancient_rank_mass_expected
        stopping_time = max((tree.time(root) for root in tree.roots), default=None)
        span = tree.span

        segments.append(
            RankSegmentSummary(
                left=tree.interval.left,
                right=tree.interval.right,
                tree_index=tree_index,
                span=span,
                root_count=root_count,
                observed_rank_count=observed_rank_count,
                missing_rank_count=missing_rank_count,
                observed_rank_mass=observed_rank_mass,
                missing_rank_mass=missing_rank_mass,
                stopping_time=stopping_time,
            )
        )

        weighted_root_count += span * root_count
        weighted_observed_rank_count += span * observed_rank_count
        weighted_observed_rank_mass += span * observed_rank_mass
        weighted_missing_rank_mass += span * missing_rank_mass
        if stopping_time is not None:
            weighted_stopping_time += span * stopping_time
            stopping_time_is_defined = True

    sequence_length = float(tree_sequence.sequence_length)
    if sequence_length <= 0.0:
        raise ValueError("Tree sequence sequence_length must be positive")

    mean_root_count = weighted_root_count / sequence_length
    mean_observed_rank_count = weighted_observed_rank_count / sequence_length
    observed_rank_mass_fraction = weighted_observed_rank_mass / (sequence_length * full_rank_mass) if full_rank_mass else 0.0
    missing_rank_mass_fraction = weighted_missing_rank_mass / (sequence_length * full_rank_mass) if full_rank_mass else 0.0
    mean_stopping_time = weighted_stopping_time / sequence_length if stopping_time_is_defined else None

    return TreeSequenceRankSummary(
        path=source_path,
        sample_count_config=sample_count_config,
        stopping_lineages_config=stopping_lineages_config,
        sample_count=sample_count,
        sequence_length=sequence_length,
        tree_count=tree_count,
        full_rank_mass=full_rank_mass,
        observed_rank_mass_expected=observed_rank_mass_expected,
        ancient_rank_mass_expected=ancient_rank_mass_expected,
        ancient_rank_value=ancient_rank_value_expected,
        segments=tuple(segments),
        mean_root_count=mean_root_count,
        mean_observed_rank_count=mean_observed_rank_count,
        mean_stopping_time=mean_stopping_time,
        observed_rank_mass_fraction=observed_rank_mass_fraction,
        missing_rank_mass_fraction=missing_rank_mass_fraction,
    )


def rank_distance(left: TreeSequenceRankSummary, right: TreeSequenceRankSummary) -> float:
    if left.sample_count != right.sample_count:
        raise ValueError(
            "Tree sequences must use the same sample count to be compared: "
            f"{left.sample_count} != {right.sample_count}"
        )
    if not math.isclose(left.sequence_length, right.sequence_length, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(
            "Tree sequences must use the same sequence length to be compared: "
            f"{left.sequence_length} != {right.sequence_length}"
        )
    if left.stopping_lineages_config != right.stopping_lineages_config:
        raise ValueError(
            "Tree sequences must use the same stopping-lineage constant to be compared: "
            f"{left.stopping_lineages_config} != {right.stopping_lineages_config}"
        )

    sequence_length = left.sequence_length
    rank_scale = left.ancient_rank_value / left.full_rank_mass if left.full_rank_mass else 0.0

    breakpoints = {0.0, sequence_length}
    breakpoints.update(segment.left for segment in left.segments)
    breakpoints.update(segment.right for segment in left.segments)
    breakpoints.update(segment.left for segment in right.segments)
    breakpoints.update(segment.right for segment in right.segments)
    ordered_breakpoints = sorted(breakpoints)

    total = 0.0
    for start, stop in zip(ordered_breakpoints, ordered_breakpoints[1:]):
        if stop <= start:
            continue
        midpoint = (start + stop) / 2.0
        left_segment = _segment_for_position(left, midpoint)
        right_segment = _segment_for_position(right, midpoint)

        feature_distance = (
            abs(left_segment.root_count - right_segment.root_count) / max(left.stopping_lineages_config, 1)
            + (0.0 if left_segment.stopping_time is None or right_segment.stopping_time is None else abs(left_segment.stopping_time - right_segment.stopping_time) * rank_scale)
        ) / 2.0
        total += feature_distance * (stop - start)

    return total / sequence_length


def analyze_tree_sequence_folder(
    folder: Path | str,
    recursive: bool = False,
    suffixes: tuple[str, ...] = (".trees",),
    sample_count_config: int | None = None,
    stopping_lineages_config: int | None = None,
) -> RankDistanceAnalysis:
    paths = discover_tree_sequence_files(folder, recursive=recursive, suffixes=suffixes)
    summaries = tuple(
        summarize_tree_sequence(
            path,
            sample_count_config=sample_count_config,
            stopping_lineages_config=stopping_lineages_config,
        )
        for path in paths
    )
    matrix = tuple(
        tuple(0.0 if i == j else rank_distance(left, right) for j, right in enumerate(summaries))
        for i, left in enumerate(summaries)
    )
    return RankDistanceAnalysis(paths=paths, summaries=summaries, matrix=matrix)


def _segment_for_position(summary: TreeSequenceRankSummary, position: float) -> RankSegmentSummary:
    for segment in summary.segments:
        if segment.left <= position < segment.right or math.isclose(position, segment.right, abs_tol=1e-12):
            return segment
    raise ValueError(f"Could not locate a tree segment at position {position}")