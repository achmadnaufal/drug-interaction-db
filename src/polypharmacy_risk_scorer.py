"""
Polypharmacy risk scorer for drug-drug interaction burden assessment.

Given a patient's medication list, this module computes a quantitative DDI
burden score by enumerating all pairwise interactions, applying severity
weights, and flagging regimens that exceed clinical polypharmacy thresholds.

Author: github.com/achmadnaufal
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Dict, Final, FrozenSet, List, Optional, Tuple

import pandas as pd

from src.main import SEVERITY_ORDER, DrugInteractionDB, _validate_drug_name


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEVERITY_WEIGHTS: Final[Dict[str, float]] = {
    "minor": 1.0,
    "moderate": 3.0,
    "major": 7.0,
    "contraindicated": 15.0,
}

#: Score at or above which a regimen is considered high-burden.
HIGH_BURDEN_THRESHOLD: Final[float] = 10.0

#: Score at or above which a regimen is considered critical.
CRITICAL_BURDEN_THRESHOLD: Final[float] = 25.0

#: Drug count at or above which polypharmacy is flagged.
POLYPHARMACY_DRUG_COUNT: Final[int] = 5


# ---------------------------------------------------------------------------
# Result dataclasses (immutable by convention — use frozen=True)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PairScore:
    """Score contribution from a single drug pair.

    Attributes:
        drug_a: Normalised name of the first drug.
        drug_b: Normalised name of the second drug.
        severity: Highest severity level found for this pair.
        weight: Numeric weight assigned to *severity*.
        interaction_count: Number of distinct interaction rows found.
        pair_score: Total score contribution (``weight × interaction_count``).
    """

    drug_a: str
    drug_b: str
    severity: str
    weight: float
    interaction_count: int
    pair_score: float


@dataclass(frozen=True)
class RiskReport:
    """Full polypharmacy risk report for a medication list.

    Attributes:
        drugs: Tuple of normalised drug names that were evaluated.
        pair_scores: Tuple of :class:`PairScore` objects, one per interacting pair.
        total_score: Sum of all pair scores.
        interacting_pairs: Number of pairs with at least one interaction.
        total_pairs: Total number of unique pairs evaluated.
        has_contraindicated: True when any pair is contraindicated.
        polypharmacy_flagged: True when the drug count meets the threshold.
        risk_level: One of ``"low"``, ``"high"``, or ``"critical"``.
    """

    drugs: Tuple[str, ...]
    pair_scores: Tuple[PairScore, ...]
    total_score: float
    interacting_pairs: int
    total_pairs: int
    has_contraindicated: bool
    polypharmacy_flagged: bool
    risk_level: str


# ---------------------------------------------------------------------------
# Core scoring function
# ---------------------------------------------------------------------------

def _highest_severity(severities: pd.Series) -> str:
    """Return the most severe level present in *severities*.

    Args:
        severities: Series of severity strings (may be mixed case).

    Returns:
        The highest severity string from :data:`~src.main.SEVERITY_ORDER`.
        Returns ``"minor"`` if no recognised level is found.
    """
    known: FrozenSet[str] = frozenset(SEVERITY_ORDER)
    normed = severities.str.lower().str.strip()
    found = [s for s in SEVERITY_ORDER if s in normed.values and s in known]
    return found[-1] if found else "minor"


def score_regimen(
    drugs: List[str],
    db: DrugInteractionDB,
    severity_weights: Optional[Dict[str, float]] = None,
    high_threshold: float = HIGH_BURDEN_THRESHOLD,
    critical_threshold: float = CRITICAL_BURDEN_THRESHOLD,
    polypharmacy_count: int = POLYPHARMACY_DRUG_COUNT,
) -> RiskReport:
    """Score a patient medication list by total DDI burden.

    Enumerates every unique pair from *drugs*, looks up each pair in *db*,
    applies severity weights, and aggregates the results into a
    :class:`RiskReport`.  All intermediate objects are new — no existing data
    is mutated.

    Args:
        drugs: List of drug names in the patient's current regimen.  Must
               contain at least one entry; each entry must be a non-empty
               string.  Single-drug lists return a zero-score report.
        db: A :class:`~src.main.DrugInteractionDB` instance with data
            already loaded via ``load_data()``.
        severity_weights: Optional mapping of severity level to numeric weight.
                          Defaults to :data:`SEVERITY_WEIGHTS`.
        high_threshold: Score at or above which risk level is ``"high"``.
                        Defaults to :data:`HIGH_BURDEN_THRESHOLD`.
        critical_threshold: Score at or above which risk level is ``"critical"``.
                            Defaults to :data:`CRITICAL_BURDEN_THRESHOLD`.
        polypharmacy_count: Minimum drug count to flag polypharmacy.
                            Defaults to :data:`POLYPHARMACY_DRUG_COUNT`.

    Returns:
        A frozen :class:`RiskReport` summarising the DDI burden for the
        given regimen.

    Raises:
        TypeError: If *drugs* is not a list, any drug name is not a string,
                   or *db* is not a :class:`~src.main.DrugInteractionDB`.
        ValueError: If *drugs* is empty or any drug name is blank.
        RuntimeError: If *db* has no data loaded.

    Example:
        >>> from src.main import DrugInteractionDB
        >>> from src.polypharmacy_risk_scorer import score_regimen
        >>> db = DrugInteractionDB()
        >>> db.load_data("demo/sample_data.csv")
        >>> report = score_regimen(["warfarin", "aspirin", "omeprazole"], db)
        >>> print(report.total_score, report.risk_level)
        7.0 low
    """
    if not isinstance(drugs, list):
        raise TypeError(f"'drugs' must be a list, got {type(drugs).__name__!r}.")
    if not drugs:
        raise ValueError("'drugs' must not be empty.")
    if not isinstance(db, DrugInteractionDB):
        raise TypeError(
            f"'db' must be a DrugInteractionDB instance, got {type(db).__name__!r}."
        )

    weights: Dict[str, float] = severity_weights if severity_weights is not None else SEVERITY_WEIGHTS

    # Validate and normalise all drug names up front (immutable list)
    normalised: Tuple[str, ...] = tuple(
        _validate_drug_name(d, f"drugs[{i}]") for i, d in enumerate(drugs)
    )

    # Deduplicate while preserving order
    seen: dict = {}
    unique_drugs: Tuple[str, ...] = tuple(
        seen.setdefault(d, d) for d in normalised if d not in seen
    )

    # Short-circuit for single drug
    if len(unique_drugs) < 2:
        return RiskReport(
            drugs=unique_drugs,
            pair_scores=(),
            total_score=0.0,
            interacting_pairs=0,
            total_pairs=0,
            has_contraindicated=False,
            polypharmacy_flagged=len(unique_drugs) >= polypharmacy_count,
            risk_level="low",
        )

    pair_scores: List[PairScore] = []
    total_pairs: int = 0

    for drug_a, drug_b in itertools.combinations(unique_drugs, 2):
        total_pairs += 1
        hits: pd.DataFrame = db.lookup(drug_a, drug_b)
        if hits.empty:
            continue

        highest: str = _highest_severity(hits["severity"])
        weight: float = weights.get(highest, 1.0)
        interaction_count: int = len(hits)
        pair_score_val: float = weight * interaction_count

        pair_scores.append(
            PairScore(
                drug_a=drug_a,
                drug_b=drug_b,
                severity=highest,
                weight=weight,
                interaction_count=interaction_count,
                pair_score=pair_score_val,
            )
        )

    total_score: float = sum(ps.pair_score for ps in pair_scores)
    has_contraindicated: bool = any(ps.severity == "contraindicated" for ps in pair_scores)
    polypharmacy_flagged: bool = len(unique_drugs) >= polypharmacy_count

    if total_score >= critical_threshold or has_contraindicated:
        risk_level = "critical"
    elif total_score >= high_threshold:
        risk_level = "high"
    else:
        risk_level = "low"

    return RiskReport(
        drugs=unique_drugs,
        pair_scores=tuple(pair_scores),
        total_score=total_score,
        interacting_pairs=len(pair_scores),
        total_pairs=total_pairs,
        has_contraindicated=has_contraindicated,
        polypharmacy_flagged=polypharmacy_flagged,
        risk_level=risk_level,
    )


def report_to_dataframe(report: RiskReport) -> pd.DataFrame:
    """Convert a :class:`RiskReport` into a tabular DataFrame of pair scores.

    Args:
        report: A :class:`RiskReport` as returned by :func:`score_regimen`.

    Returns:
        DataFrame with columns ``drug_a``, ``drug_b``, ``severity``,
        ``weight``, ``interaction_count``, and ``pair_score``, sorted by
        ``pair_score`` descending.  Returns an empty DataFrame when no
        interacting pairs were found.

    Raises:
        TypeError: If *report* is not a :class:`RiskReport`.

    Example:
        >>> df = report_to_dataframe(report)
        >>> print(df[["drug_a", "drug_b", "severity", "pair_score"]])
    """
    if not isinstance(report, RiskReport):
        raise TypeError(
            f"'report' must be a RiskReport instance, got {type(report).__name__!r}."
        )
    if not report.pair_scores:
        return pd.DataFrame(
            columns=["drug_a", "drug_b", "severity", "weight", "interaction_count", "pair_score"]
        )

    rows = [
        {
            "drug_a": ps.drug_a,
            "drug_b": ps.drug_b,
            "severity": ps.severity,
            "weight": ps.weight,
            "interaction_count": ps.interaction_count,
            "pair_score": ps.pair_score,
        }
        for ps in report.pair_scores
    ]
    df = pd.DataFrame(rows)
    return df.sort_values("pair_score", ascending=False).reset_index(drop=True)
