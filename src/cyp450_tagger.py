"""
CYP450 enzyme-mediated interaction tagger.

This module scans the ``mechanism`` text of drug-drug interaction records and
tags each row with the cytochrome P450 isoenzymes (and related transporters
such as P-glycoprotein) implicated in the interaction.  The output is useful
for filtering pharmacokinetic interactions by metabolic pathway and for
generating clinician-facing summaries grouped by enzyme.

Note:
    This module is for educational and illustrative use only.  It does not
    constitute clinical advice and is not endorsed by the FDA or any other
    regulatory body.

Author: github.com/achmadnaufal
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Final, FrozenSet, List, Optional, Tuple

import pandas as pd

from src.main import DrugInteractionDB


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Recognised CYP450 isoenzymes and related drug-handling proteins.  Order is
#: significant for deterministic tagging output.
KNOWN_ENZYMES: Final[Tuple[str, ...]] = (
    "CYP1A2",
    "CYP2B6",
    "CYP2C8",
    "CYP2C9",
    "CYP2C19",
    "CYP2D6",
    "CYP2E1",
    "CYP3A4",
    "CYP3A5",
    "P-glycoprotein",
    "UGT1A1",
    "OATP1B1",
)

#: Words in the ``mechanism`` field that signal an inhibitor relationship.
INHIBITOR_KEYWORDS: Final[FrozenSet[str]] = frozenset(
    {"inhibit", "inhibits", "inhibition", "inhibitor", "inhibiting"}
)

#: Words in the ``mechanism`` field that signal an inducer relationship.
INDUCER_KEYWORDS: Final[FrozenSet[str]] = frozenset(
    {"induce", "induces", "induction", "inducer", "inducing"}
)

#: Sentinel role string when neither inhibition nor induction is detected.
ROLE_UNKNOWN: Final[str] = "unspecified"


# ---------------------------------------------------------------------------
# Result dataclass (immutable)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EnzymeTag:
    """Tag describing one enzyme involvement in a single interaction row.

    Attributes:
        enzyme: Canonical enzyme identifier (e.g. ``"CYP3A4"``).
        role: One of ``"inhibitor"``, ``"inducer"``, or ``"unspecified"``.
    """

    enzyme: str
    role: str


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def _detect_role(mechanism: str) -> str:
    """Detect the metabolic role described in *mechanism* text.

    Args:
        mechanism: Free-text mechanism description (case-insensitive).

    Returns:
        ``"inhibitor"`` if any inhibitor keyword is present, ``"inducer"``
        if any inducer keyword is present, otherwise :data:`ROLE_UNKNOWN`.
        Inhibition takes precedence when both are mentioned.
    """
    text = mechanism.lower()
    tokens = set(re.findall(r"[a-z]+", text))
    if tokens & INHIBITOR_KEYWORDS:
        return "inhibitor"
    if tokens & INDUCER_KEYWORDS:
        return "inducer"
    return ROLE_UNKNOWN


def extract_enzymes(
    mechanism: str,
    enzymes: Optional[Tuple[str, ...]] = None,
) -> Tuple[EnzymeTag, ...]:
    """Extract enzyme tags from a single mechanism description.

    Matching is case-insensitive and word-boundary aware so that ``CYP3A4``
    is not falsely matched inside an unrelated token.  Each enzyme appears at
    most once in the returned tuple regardless of how many times it occurs in
    the mechanism text.

    Args:
        mechanism: Free-text mechanism description from a DDI record.  Empty
                   or non-string input yields an empty tuple.
        enzymes: Optional iterable of enzyme identifiers to look for.
                 Defaults to :data:`KNOWN_ENZYMES`.

    Returns:
        Tuple of :class:`EnzymeTag` objects in the order defined by *enzymes*.

    Raises:
        TypeError: If *mechanism* is not a string.
    """
    if mechanism is None or (isinstance(mechanism, float) and pd.isna(mechanism)):
        return ()
    if not isinstance(mechanism, str):
        raise TypeError(
            f"'mechanism' must be a string, got {type(mechanism).__name__!r}."
        )

    if not mechanism.strip():
        return ()

    catalogue: Tuple[str, ...] = enzymes if enzymes is not None else KNOWN_ENZYMES
    role: str = _detect_role(mechanism)
    text = mechanism.lower()

    tags: List[EnzymeTag] = []
    seen: set = set()
    for enzyme in catalogue:
        pattern = r"(?<![a-z0-9])" + re.escape(enzyme.lower()) + r"(?![a-z0-9])"
        if re.search(pattern, text) and enzyme not in seen:
            tags.append(EnzymeTag(enzyme=enzyme, role=role))
            seen.add(enzyme)
    return tuple(tags)


def tag_interactions(
    df: pd.DataFrame,
    enzymes: Optional[Tuple[str, ...]] = None,
) -> pd.DataFrame:
    """Append CYP450 enzyme tags to a DDI DataFrame.

    The input DataFrame is never mutated; a new DataFrame is returned with
    two additional columns:

    - ``cyp_enzymes``: comma-separated list of enzymes mentioned, or empty.
    - ``cyp_role``: ``"inhibitor"``, ``"inducer"``, ``"unspecified"``, or
                    empty string when no enzyme is detected.

    Args:
        df: DataFrame containing at least a ``mechanism`` column.
        enzymes: Optional enzyme catalogue.  Defaults to :data:`KNOWN_ENZYMES`.

    Returns:
        New DataFrame with the two extra columns appended.

    Raises:
        TypeError: If *df* is not a pandas DataFrame.
        ValueError: If the ``mechanism`` column is missing.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"'df' must be a pandas DataFrame, got {type(df).__name__!r}."
        )
    if "mechanism" not in df.columns:
        raise ValueError(
            "Input DataFrame must contain a 'mechanism' column for CYP tagging."
        )

    tag_lists: List[Tuple[EnzymeTag, ...]] = [
        extract_enzymes(m, enzymes=enzymes) for m in df["mechanism"].tolist()
    ]
    enzyme_strs: List[str] = [
        ", ".join(t.enzyme for t in tags) for tags in tag_lists
    ]
    role_strs: List[str] = [
        tags[0].role if tags else "" for tags in tag_lists
    ]

    return df.assign(cyp_enzymes=enzyme_strs, cyp_role=role_strs)


def summarise_by_enzyme(
    df: pd.DataFrame,
    enzymes: Optional[Tuple[str, ...]] = None,
) -> pd.DataFrame:
    """Summarise interaction counts grouped by CYP enzyme.

    Args:
        df: DataFrame containing a ``mechanism`` column (raw or pre-tagged).
        enzymes: Optional enzyme catalogue.  Defaults to :data:`KNOWN_ENZYMES`.

    Returns:
        DataFrame with columns ``enzyme``, ``inhibitor_count``,
        ``inducer_count``, ``unspecified_count``, and ``total`` sorted by
        ``total`` descending.  Returns an empty DataFrame with the expected
        columns when no enzyme mentions are present.

    Raises:
        TypeError: If *df* is not a pandas DataFrame.
        ValueError: If the ``mechanism`` column is missing.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"'df' must be a pandas DataFrame, got {type(df).__name__!r}."
        )
    if "mechanism" not in df.columns:
        raise ValueError(
            "Input DataFrame must contain a 'mechanism' column for summarisation."
        )

    counts: Dict[str, Dict[str, int]] = {}
    for mechanism in df["mechanism"].tolist():
        for tag in extract_enzymes(mechanism, enzymes=enzymes):
            bucket = counts.setdefault(
                tag.enzyme,
                {"inhibitor": 0, "inducer": 0, ROLE_UNKNOWN: 0},
            )
            bucket[tag.role] += 1

    columns = ["enzyme", "inhibitor_count", "inducer_count", "unspecified_count", "total"]
    if not counts:
        return pd.DataFrame(columns=columns)

    rows = [
        {
            "enzyme": enzyme,
            "inhibitor_count": bucket["inhibitor"],
            "inducer_count": bucket["inducer"],
            "unspecified_count": bucket[ROLE_UNKNOWN],
            "total": bucket["inhibitor"] + bucket["inducer"] + bucket[ROLE_UNKNOWN],
        }
        for enzyme, bucket in counts.items()
    ]
    out = pd.DataFrame(rows, columns=columns)
    return out.sort_values(["total", "enzyme"], ascending=[False, True]).reset_index(drop=True)


def filter_by_enzyme(
    db: DrugInteractionDB,
    enzyme: str,
    role: Optional[str] = None,
) -> pd.DataFrame:
    """Return interactions in *db* that mention *enzyme* in their mechanism.

    Args:
        db: A :class:`~src.main.DrugInteractionDB` with data loaded.
        enzyme: Enzyme identifier (case-insensitive), e.g. ``"CYP3A4"``.
        role: Optional role filter — one of ``"inhibitor"``, ``"inducer"``,
              or ``"unspecified"``.  When omitted all roles are returned.

    Returns:
        DataFrame of matching interaction rows with ``cyp_enzymes`` and
        ``cyp_role`` columns appended.

    Raises:
        TypeError: If *db* is not a :class:`DrugInteractionDB` or *enzyme*
                   is not a string.
        ValueError: If *enzyme* is empty/unknown or *role* is not a
                    recognised value.
        RuntimeError: If *db* has no data loaded.
    """
    if not isinstance(db, DrugInteractionDB):
        raise TypeError(
            f"'db' must be a DrugInteractionDB instance, got {type(db).__name__!r}."
        )
    if not isinstance(enzyme, str):
        raise TypeError(
            f"'enzyme' must be a string, got {type(enzyme).__name__!r}."
        )
    cleaned = enzyme.strip()
    if not cleaned:
        raise ValueError("'enzyme' must not be empty or whitespace-only.")

    canonical_lookup = {e.lower(): e for e in KNOWN_ENZYMES}
    canonical = canonical_lookup.get(cleaned.lower())
    if canonical is None:
        raise ValueError(
            f"Unknown enzyme '{enzyme}'. Choose from: {list(KNOWN_ENZYMES)}"
        )

    valid_roles = {"inhibitor", "inducer", ROLE_UNKNOWN}
    if role is not None and role not in valid_roles:
        raise ValueError(
            f"Unknown role '{role}'. Choose from: {sorted(valid_roles)}"
        )

    source = db._require_loaded()
    tagged = tag_interactions(source)

    mask = tagged["cyp_enzymes"].str.contains(
        rf"(?:^|, ){re.escape(canonical)}(?:,|$)", regex=True, na=False
    )
    filtered = tagged.loc[mask].copy()
    if role is not None:
        filtered = filtered.loc[filtered["cyp_role"] == role].copy()
    return filtered.reset_index(drop=True)
