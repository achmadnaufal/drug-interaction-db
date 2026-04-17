"""
Tests for src/cyp450_tagger.py.

Run with:
    pytest tests/test_cyp450_tagger.py -v
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.cyp450_tagger import (
    KNOWN_ENZYMES,
    ROLE_UNKNOWN,
    EnzymeTag,
    extract_enzymes,
    filter_by_enzyme,
    summarise_by_enzyme,
    tag_interactions,
)
from src.main import DrugInteractionDB

SAMPLE_CSV = str(Path(__file__).parent.parent / "demo" / "sample_data.csv")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def loaded_db() -> DrugInteractionDB:
    """Return a DrugInteractionDB loaded with the sample dataset."""
    db = DrugInteractionDB()
    db.load_data(SAMPLE_CSV)
    return db


@pytest.fixture
def small_df() -> pd.DataFrame:
    """Tiny DataFrame for unit-level tagging tests."""
    return pd.DataFrame(
        [
            {
                "drug_a": "warfarin",
                "drug_b": "fluconazole",
                "interaction_type": "pharmacokinetic",
                "severity": "major",
                "mechanism": "Fluconazole inhibits CYP2C9 reducing warfarin metabolism",
            },
            {
                "drug_a": "rifampicin",
                "drug_b": "oral_contraceptives",
                "interaction_type": "pharmacokinetic",
                "severity": "major",
                "mechanism": "Rifampicin induces CYP3A4 increasing hormone metabolism",
            },
            {
                "drug_a": "warfarin",
                "drug_b": "aspirin",
                "interaction_type": "pharmacodynamic",
                "severity": "major",
                "mechanism": "Additive anticoagulant and antiplatelet effects",
            },
        ]
    )


# ---------------------------------------------------------------------------
# 1. Happy path — extract_enzymes
# ---------------------------------------------------------------------------

class TestExtractEnzymes:
    def test_inhibitor_role_detected(self) -> None:
        tags = extract_enzymes("Fluconazole inhibits CYP2C9 reducing warfarin metabolism")
        assert tags == (EnzymeTag(enzyme="CYP2C9", role="inhibitor"),)

    def test_inducer_role_detected(self) -> None:
        tags = extract_enzymes("Rifampicin induces CYP3A4 increasing metabolism")
        assert tags == (EnzymeTag(enzyme="CYP3A4", role="inducer"),)

    def test_multiple_enzymes_in_one_mechanism(self) -> None:
        tags = extract_enzymes(
            "Amiodarone inhibits CYP3A4 and CYP2C9 reducing statin clearance"
        )
        enzymes = [t.enzyme for t in tags]
        assert "CYP3A4" in enzymes
        assert "CYP2C9" in enzymes
        assert all(t.role == "inhibitor" for t in tags)

    def test_pglycoprotein_detected(self) -> None:
        tags = extract_enzymes(
            "Amiodarone inhibits P-glycoprotein and renal tubular secretion"
        )
        assert any(t.enzyme == "P-glycoprotein" for t in tags)

    def test_case_insensitive_match(self) -> None:
        tags = extract_enzymes("INDUCES cyp3a4 markedly")
        assert tags == (EnzymeTag(enzyme="CYP3A4", role="inducer"),)


# ---------------------------------------------------------------------------
# 2. Edge cases — extract_enzymes
# ---------------------------------------------------------------------------

class TestExtractEnzymesEdges:
    def test_empty_string_returns_empty_tuple(self) -> None:
        assert extract_enzymes("") == ()

    def test_whitespace_string_returns_empty_tuple(self) -> None:
        assert extract_enzymes("   \t\n  ") == ()

    def test_none_returns_empty_tuple(self) -> None:
        assert extract_enzymes(None) == ()  # type: ignore[arg-type]

    def test_nan_returns_empty_tuple(self) -> None:
        import math
        assert extract_enzymes(math.nan) == ()  # type: ignore[arg-type]

    def test_no_enzymes_mentioned(self) -> None:
        assert extract_enzymes("Additive anticoagulant and antiplatelet effects") == ()

    def test_unspecified_role_when_no_keyword(self) -> None:
        tags = extract_enzymes("CYP3A4 substrate interaction noted")
        assert tags == (EnzymeTag(enzyme="CYP3A4", role=ROLE_UNKNOWN),)

    def test_inhibition_takes_precedence_over_induction(self) -> None:
        tags = extract_enzymes("CYP3A4 inducer that paradoxically inhibits at low dose")
        assert tags[0].role == "inhibitor"

    def test_word_boundary_avoids_false_match(self) -> None:
        # 'XCYP3A4Y' should not be treated as CYP3A4 due to word-boundary guard
        tags = extract_enzymes("inhibits xcyp3a4y compound")
        assert tags == ()

    def test_non_string_mechanism_raises_typeerror(self) -> None:
        with pytest.raises(TypeError):
            extract_enzymes(123)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# 3. tag_interactions DataFrame helper
# ---------------------------------------------------------------------------

class TestTagInteractions:
    def test_appends_two_columns(self, small_df: pd.DataFrame) -> None:
        out = tag_interactions(small_df)
        assert "cyp_enzymes" in out.columns
        assert "cyp_role" in out.columns

    def test_does_not_mutate_input(self, small_df: pd.DataFrame) -> None:
        before_cols = list(small_df.columns)
        _ = tag_interactions(small_df)
        assert list(small_df.columns) == before_cols

    def test_pharmacodynamic_row_has_no_enzymes(self, small_df: pd.DataFrame) -> None:
        out = tag_interactions(small_df)
        pd_row = out.loc[out["drug_a"] == "warfarin"].iloc[-1]
        assert pd_row["cyp_enzymes"] == ""
        assert pd_row["cyp_role"] == ""

    def test_inducer_row_tagged_correctly(self, small_df: pd.DataFrame) -> None:
        out = tag_interactions(small_df)
        rifampicin_row = out.loc[out["drug_a"] == "rifampicin"].iloc[0]
        assert rifampicin_row["cyp_enzymes"] == "CYP3A4"
        assert rifampicin_row["cyp_role"] == "inducer"

    def test_missing_mechanism_column_raises(self) -> None:
        df = pd.DataFrame({"drug_a": ["a"], "drug_b": ["b"]})
        with pytest.raises(ValueError, match="mechanism"):
            tag_interactions(df)

    def test_non_dataframe_raises_typeerror(self) -> None:
        with pytest.raises(TypeError):
            tag_interactions("not a dataframe")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# 4. summarise_by_enzyme aggregation
# ---------------------------------------------------------------------------

class TestSummariseByEnzyme:
    def test_summary_columns(self, loaded_db: DrugInteractionDB) -> None:
        summary = summarise_by_enzyme(loaded_db._df)
        assert list(summary.columns) == [
            "enzyme",
            "inhibitor_count",
            "inducer_count",
            "unspecified_count",
            "total",
        ]

    def test_sample_data_includes_cyp3a4(self, loaded_db: DrugInteractionDB) -> None:
        summary = summarise_by_enzyme(loaded_db._df)
        assert "CYP3A4" in summary["enzyme"].tolist()

    def test_sorted_by_total_descending(self, loaded_db: DrugInteractionDB) -> None:
        summary = summarise_by_enzyme(loaded_db._df)
        totals = summary["total"].tolist()
        assert totals == sorted(totals, reverse=True)

    def test_empty_dataframe_returns_empty_summary(self) -> None:
        df = pd.DataFrame(columns=["mechanism"])
        summary = summarise_by_enzyme(df)
        assert summary.empty
        assert "enzyme" in summary.columns

    def test_no_enzyme_mentions_returns_empty(self) -> None:
        df = pd.DataFrame({"mechanism": ["additive effects", "synergy noted"]})
        summary = summarise_by_enzyme(df)
        assert summary.empty


# ---------------------------------------------------------------------------
# 5. filter_by_enzyme integration
# ---------------------------------------------------------------------------

class TestFilterByEnzyme:
    def test_filter_known_enzyme_returns_rows(self, loaded_db: DrugInteractionDB) -> None:
        result = filter_by_enzyme(loaded_db, "CYP3A4")
        assert not result.empty
        assert all("CYP3A4" in tags for tags in result["cyp_enzymes"])

    def test_filter_case_insensitive(self, loaded_db: DrugInteractionDB) -> None:
        upper = filter_by_enzyme(loaded_db, "CYP3A4")
        lower = filter_by_enzyme(loaded_db, "cyp3a4")
        assert len(upper) == len(lower)

    def test_role_filter_inhibitor_only(self, loaded_db: DrugInteractionDB) -> None:
        result = filter_by_enzyme(loaded_db, "CYP3A4", role="inhibitor")
        assert all(r == "inhibitor" for r in result["cyp_role"])

    def test_unknown_enzyme_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError, match="Unknown enzyme"):
            filter_by_enzyme(loaded_db, "CYP9Z9")

    def test_empty_enzyme_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError, match="empty"):
            filter_by_enzyme(loaded_db, "   ")

    def test_invalid_role_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError, match="Unknown role"):
            filter_by_enzyme(loaded_db, "CYP3A4", role="bogus")

    def test_non_db_raises_typeerror(self) -> None:
        with pytest.raises(TypeError):
            filter_by_enzyme("not a db", "CYP3A4")  # type: ignore[arg-type]

    def test_non_string_enzyme_raises_typeerror(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(TypeError):
            filter_by_enzyme(loaded_db, 123)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# 6. Constants sanity
# ---------------------------------------------------------------------------

class TestConstants:
    def test_known_enzymes_are_unique(self) -> None:
        assert len(KNOWN_ENZYMES) == len(set(KNOWN_ENZYMES))

    def test_known_enzymes_contains_core_cyps(self) -> None:
        for required in ("CYP3A4", "CYP2C9", "CYP2C19", "CYP2D6", "CYP1A2"):
            assert required in KNOWN_ENZYMES

    def test_p_glycoprotein_present(self) -> None:
        assert "P-glycoprotein" in KNOWN_ENZYMES
