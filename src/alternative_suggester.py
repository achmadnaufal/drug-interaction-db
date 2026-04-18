"""
Alternative-drug suggester for drug-drug interaction mitigation.

Given a target drug that is causing a problematic interaction with a
patient's existing regimen, this module proposes candidate replacement
drugs from the same therapeutic class that carry a lower overall DDI
burden against the rest of the regimen.

Core ideas:

1. A curated, immutable ``THERAPEUTIC_CLASSES`` mapping groups drugs with
   similar clinical indications (e.g. PPIs, SSRIs, statins).
2. :func:`suggest_alternatives` scores each candidate in the same class
   against every other drug in the patient's regimen using the severity
   weights from the polypharmacy scorer.
3. Results are returned as immutable :class:`AlternativeSuggestion`
   tuples, sorted by ascending DDI burden (the safest alternatives first).

All functions are pure — neither the supplied DataFrame nor the
supplied regimen list is mutated.  Unknown drugs, empty regimens,
duplicate entries, and case variations are all handled explicitly.

Author: github.com/achmadnaufal
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Final, FrozenSet, List, Mapping, Optional, Tuple

import pandas as pd

from src.main import DrugInteractionDB, _validate_drug_name
from src.polypharmacy_risk_scorer import SEVERITY_WEIGHTS, _highest_severity


# ---------------------------------------------------------------------------
# Curated therapeutic class registry (immutable)
# ---------------------------------------------------------------------------

#: Mapping of therapeutic class name (lowercase) to the tuple of drugs
#: that belong to that class.  Entries are deliberately kept lowercase
#: so lookups stay case-insensitive without extra normalisation.
THERAPEUTIC_CLASSES: Final[Mapping[str, Tuple[str, ...]]] = {
    "ppi": (
        "omeprazole",
        "esomeprazole",
        "pantoprazole",
        "lansoprazole",
        "rabeprazole",
    ),
    "ssri": (
        "fluoxetine",
        "sertraline",
        "citalopram",
        "escitalopram",
        "paroxetine",
    ),
    "statin": (
        "simvastatin",
        "atorvastatin",
        "rosuvastatin",
        "pravastatin",
        "fluvastatin",
        "lovastatin",
    ),
    "macrolide": (
        "clarithromycin",
        "erythromycin",
        "azithromycin",
    ),
    "azole_antifungal": (
        "fluconazole",
        "itraconazole",
        "ketoconazole",
        "voriconazole",
        "posaconazole",
    ),
    "nsaid": (
        "ibuprofen",
        "naproxen",
        "diclofenac",
        "celecoxib",
        "ketorolac",
    ),
    "ace_inhibitor": (
        "lisinopril",
        "enalapril",
        "ramipril",
        "captopril",
        "benazepril",
    ),
    "anticoagulant": (
        "warfarin",
        "apixaban",
        "rivaroxaban",
        "dabigatran",
        "edoxaban",
    ),
}


# ---------------------------------------------------------------------------
# Result dataclass (frozen / immutable)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AlternativeSuggestion:
    """A single candidate alternative and its DDI burden versus the regimen.

    Attributes:
        candidate: Normalised name of the proposed alternative drug.
        therapeutic_class: Class name the candidate shares with the target.
        total_score: Sum of severity weights for all interactions between
            *candidate* and every other drug in the evaluated regimen.
        highest_severity: Most severe interaction level found, or
            ``"none"`` when the candidate has no interactions with the
            regimen.
        interacting_drugs: Tuple of regimen drugs (other than the
            replaced target) that interact with *candidate*.
        is_safer: ``True`` when ``total_score`` is strictly lower than the
            baseline score of the drug being replaced.
    """

    candidate: str
    therapeutic_class: str
    total_score: float
    highest_severity: str
    interacting_drugs: Tuple[str, ...]
    is_safer: bool


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _find_class(drug: str) -> Optional[str]:
    """Return the therapeutic class containing *drug*, or ``None``.

    Matching is case-insensitive and trims whitespace.  When a drug
    belongs to multiple classes (unusual but possible in pharmacology),
    the first registered class wins — the registry is ordered.

    Args:
        drug: Normalised drug name.

    Returns:
        The lowercase class name, or ``None`` when *drug* is not a
        member of any registered therapeutic class.
    """
    normalised = drug.strip().lower()
    for class_name, members in THERAPEUTIC_CLASSES.items():
        if normalised in members:
            return class_name
    return None


def _score_candidate(
    candidate: str,
    regimen: Tuple[str, ...],
    db: DrugInteractionDB,
    weights: Mapping[str, float],
) -> Tuple[float, str, Tuple[str, ...]]:
    """Compute the DDI burden for *candidate* against *regimen*.

    Args:
        candidate: Normalised candidate drug name.
        regimen: Tuple of other regimen drugs to evaluate against.
        db: Loaded :class:`DrugInteractionDB` instance.
        weights: Severity-to-weight mapping.

    Returns:
        Tuple of ``(total_score, highest_severity, interacting_drugs)``.
        ``highest_severity`` is ``"none"`` when the candidate has zero
        recorded interactions with the regimen.
    """
    score: float = 0.0
    highest: str = "none"
    interacting: List[str] = []

    for other in regimen:
        if other == candidate:
            continue
        hits: pd.DataFrame = db.lookup(candidate, other)
        if hits.empty:
            continue

        pair_sev = _highest_severity(hits["severity"])
        score += weights.get(pair_sev, 0.0) * len(hits)
        interacting.append(other)

        # Track the overall worst severity across all pairs.
        if highest == "none":
            highest = pair_sev
        else:
            # Pick the higher-ordinal severity between `highest` and pair_sev.
            order = ["minor", "moderate", "major", "contraindicated"]
            if order.index(pair_sev) > order.index(highest):
                highest = pair_sev

    return score, highest, tuple(interacting)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def suggest_alternatives(
    target: str,
    regimen: List[str],
    db: DrugInteractionDB,
    severity_weights: Optional[Dict[str, float]] = None,
    max_results: int = 5,
) -> Tuple[AlternativeSuggestion, ...]:
    """Suggest safer alternatives for *target* given the patient *regimen*.

    The function identifies the therapeutic class of *target*, scores every
    other class member against the remainder of the regimen, and returns
    the candidates with the lowest total DDI burden.

    Args:
        target: The drug to replace.  Must appear in one of the
            :data:`THERAPEUTIC_CLASSES` groups; otherwise an empty tuple
            is returned.  Matching is case-insensitive.
        regimen: Full current medication list, including *target*.  Must
            contain at least one entry.  Each entry must be a non-empty
            string.  The list is deduplicated internally.
        db: A loaded :class:`DrugInteractionDB` instance.
        severity_weights: Optional override for severity-to-weight mapping.
            Defaults to :data:`~src.polypharmacy_risk_scorer.SEVERITY_WEIGHTS`.
        max_results: Maximum number of suggestions to return (top-N by
            ascending total score).  Must be a positive integer.

    Returns:
        Tuple of :class:`AlternativeSuggestion` objects, sorted by
        ascending ``total_score`` (safest first).  Returns an empty tuple
        when *target* is not in any registered class, when the class has
        no other members, or when no candidates exist.

    Raises:
        TypeError: If *target* is not a string, *regimen* is not a list,
            *db* is not a :class:`DrugInteractionDB`, or *max_results*
            is not an integer.
        ValueError: If *target* is empty, *regimen* is empty, any entry
            in *regimen* is blank, or *max_results* is non-positive.
        RuntimeError: If *db* has no data loaded.

    Example:
        >>> from src.main import DrugInteractionDB
        >>> from src.alternative_suggester import suggest_alternatives
        >>> db = DrugInteractionDB()
        >>> db.load_data("demo/sample_data.csv")
        >>> suggestions = suggest_alternatives(
        ...     target="omeprazole",
        ...     regimen=["clopidogrel", "omeprazole", "aspirin"],
        ...     db=db,
        ... )
        >>> for s in suggestions:
        ...     print(s.candidate, s.total_score, s.highest_severity)
    """
    # --- Input validation (fail fast, explicit messages) ---
    target_validated = _validate_drug_name(target, "target")

    if not isinstance(regimen, list):
        raise TypeError(
            f"'regimen' must be a list, got {type(regimen).__name__!r}."
        )
    if not regimen:
        raise ValueError("'regimen' must not be empty.")
    if not isinstance(db, DrugInteractionDB):
        raise TypeError(
            f"'db' must be a DrugInteractionDB instance, got {type(db).__name__!r}."
        )
    if not isinstance(max_results, int) or isinstance(max_results, bool):
        raise TypeError(
            f"'max_results' must be an int, got {type(max_results).__name__!r}."
        )
    if max_results <= 0:
        raise ValueError(
            f"'max_results' must be positive, got {max_results}."
        )

    weights: Mapping[str, float] = (
        severity_weights if severity_weights is not None else SEVERITY_WEIGHTS
    )

    # Normalise and deduplicate the regimen (preserves first-seen order).
    normalised: List[str] = [
        _validate_drug_name(d, f"regimen[{i}]") for i, d in enumerate(regimen)
    ]
    seen: FrozenSet[str] = frozenset()
    unique_regimen: List[str] = []
    for d in normalised:
        if d not in seen:
            unique_regimen.append(d)
            seen = seen | {d}

    # Identify target's therapeutic class; abort cleanly if unknown.
    class_name = _find_class(target_validated)
    if class_name is None:
        return ()

    # Build the "rest of regimen" — everything except the target.
    rest: Tuple[str, ...] = tuple(d for d in unique_regimen if d != target_validated)

    # Score the target itself so we can flag safer alternatives.
    baseline_score, _, _ = _score_candidate(
        target_validated, rest, db, weights
    )

    # Score each same-class candidate.
    candidates: Tuple[str, ...] = THERAPEUTIC_CLASSES[class_name]
    suggestions: List[AlternativeSuggestion] = []

    for cand in candidates:
        if cand == target_validated:
            continue
        score, highest, interacting = _score_candidate(cand, rest, db, weights)
        suggestions.append(
            AlternativeSuggestion(
                candidate=cand,
                therapeutic_class=class_name,
                total_score=score,
                highest_severity=highest,
                interacting_drugs=interacting,
                is_safer=score < baseline_score,
            )
        )

    # Sort by ascending score (safest first), then alphabetical for
    # deterministic output under ties.
    suggestions.sort(key=lambda s: (s.total_score, s.candidate))

    return tuple(suggestions[:max_results])


def suggestions_to_dataframe(
    suggestions: Tuple[AlternativeSuggestion, ...],
) -> pd.DataFrame:
    """Flatten a tuple of suggestions into a tabular DataFrame.

    Args:
        suggestions: Tuple of :class:`AlternativeSuggestion` objects as
            returned by :func:`suggest_alternatives`.

    Returns:
        DataFrame with columns ``candidate``, ``therapeutic_class``,
        ``total_score``, ``highest_severity``, ``interacting_drugs``,
        and ``is_safer``.  Returns an empty DataFrame when
        *suggestions* is empty.

    Raises:
        TypeError: If *suggestions* is not a tuple or contains non-
            :class:`AlternativeSuggestion` items.
    """
    if not isinstance(suggestions, tuple):
        raise TypeError(
            f"'suggestions' must be a tuple, got {type(suggestions).__name__!r}."
        )

    columns = [
        "candidate",
        "therapeutic_class",
        "total_score",
        "highest_severity",
        "interacting_drugs",
        "is_safer",
    ]
    if not suggestions:
        return pd.DataFrame(columns=columns)

    rows = []
    for s in suggestions:
        if not isinstance(s, AlternativeSuggestion):
            raise TypeError(
                "Every item in 'suggestions' must be an AlternativeSuggestion, "
                f"got {type(s).__name__!r}."
            )
        rows.append(
            {
                "candidate": s.candidate,
                "therapeutic_class": s.therapeutic_class,
                "total_score": s.total_score,
                "highest_severity": s.highest_severity,
                "interacting_drugs": ", ".join(s.interacting_drugs),
                "is_safer": s.is_safer,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def list_classes() -> Tuple[str, ...]:
    """Return the tuple of therapeutic class names currently registered.

    Returns:
        Tuple of lowercase class names in registration order.
    """
    return tuple(THERAPEUTIC_CLASSES.keys())
