"""
Unit tests for DrugInteractionDB.

Run with:
    pytest tests/test_interactions.py -v

Coverage targets:
  - Interaction lookup (exact and bidirectional)
  - Severity filtering (exact and minimum-threshold)
  - Bidirectional search (drug_a <-> drug_b symmetry)
  - Multi-drug interaction check
  - Edge cases: unknown drug, self-interaction, empty result
  - Case-insensitive matching
  - Input validation and error handling
"""

import pytest
import pandas as pd
import sys
from pathlib import Path

# Ensure project root is importable regardless of where pytest is invoked
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.main import DrugInteractionDB, SEVERITY_ORDER


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_ROWS = [
    {
        "drug_a": "warfarin",
        "drug_b": "aspirin",
        "interaction_type": "pharmacodynamic",
        "severity": "major",
        "mechanism": "Additive anticoagulant effects",
        "clinical_effect": "Increased bleeding risk",
        "recommendation": "Avoid combination",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "simvastatin",
        "drug_b": "amiodarone",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "CYP3A4 inhibition",
        "clinical_effect": "Myopathy risk",
        "recommendation": "Limit simvastatin dose",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "clopidogrel",
        "drug_b": "omeprazole",
        "interaction_type": "pharmacokinetic",
        "severity": "moderate",
        "mechanism": "CYP2C19 inhibition",
        "clinical_effect": "Reduced antiplatelet effect",
        "recommendation": "Use pantoprazole",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "ssris",
        "drug_b": "maois",
        "interaction_type": "pharmacodynamic",
        "severity": "contraindicated",
        "mechanism": "Serotonin syndrome",
        "clinical_effect": "Life-threatening serotonin syndrome",
        "recommendation": "Absolutely contraindicated",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "carbamazepine",
        "drug_b": "valproate",
        "interaction_type": "pharmacokinetic",
        "severity": "moderate",
        "mechanism": "Mutual enzyme induction",
        "clinical_effect": "Unpredictable drug levels",
        "recommendation": "Monitor serum levels",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "lithium",
        "drug_b": "ibuprofen",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Reduced renal clearance",
        "clinical_effect": "Lithium toxicity",
        "recommendation": "Avoid NSAIDs",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "sildenafil",
        "drug_b": "nitrates",
        "interaction_type": "pharmacodynamic",
        "severity": "contraindicated",
        "mechanism": "Additive vasodilation",
        "clinical_effect": "Severe hypotension",
        "recommendation": "Absolute contraindication",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "digoxin",
        "drug_b": "amiodarone",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "P-glycoprotein inhibition",
        "clinical_effect": "Digoxin toxicity",
        "recommendation": "Reduce digoxin 50%",
        "evidence_level": "A",
        "source": "Test",
    },
    {
        "drug_a": "metformin",
        "drug_b": "contrast_iodine",
        "interaction_type": "pharmacodynamic",
        "severity": "major",
        "mechanism": "Nephrotoxicity-mediated accumulation",
        "clinical_effect": "Lactic acidosis risk",
        "recommendation": "Hold metformin 48 h",
        "evidence_level": "B",
        "source": "Test",
    },
    {
        "drug_a": "phenytoin",
        "drug_b": "folic_acid",
        "interaction_type": "pharmacokinetic",
        "severity": "minor",
        "mechanism": "Increased phenytoin metabolism",
        "clinical_effect": "Reduced phenytoin levels",
        "recommendation": "Monitor phenytoin levels",
        "evidence_level": "B",
        "source": "Test",
    },
]


@pytest.fixture()
def sample_df() -> pd.DataFrame:
    """Return a small interaction DataFrame for testing."""
    return pd.DataFrame(SAMPLE_ROWS)


@pytest.fixture()
def db(sample_df: pd.DataFrame) -> DrugInteractionDB:
    """Return a DrugInteractionDB instance pre-loaded with sample data."""
    instance = DrugInteractionDB()
    instance._df = instance.preprocess(sample_df)
    return instance


@pytest.fixture()
def db_empty() -> DrugInteractionDB:
    """Return a DrugInteractionDB instance with no data loaded."""
    return DrugInteractionDB()


# ---------------------------------------------------------------------------
# Test 1: Interaction lookup — known pair
# ---------------------------------------------------------------------------

class TestLookup:
    def test_lookup_known_pair_returns_rows(self, db: DrugInteractionDB) -> None:
        """lookup() returns at least one row for a known drug pair."""
        result = db.lookup("warfarin", "aspirin")
        assert not result.empty, "Expected interaction rows for warfarin/aspirin"

    def test_lookup_known_pair_correct_severity(self, db: DrugInteractionDB) -> None:
        """Returned row for warfarin/aspirin has severity 'major'."""
        result = db.lookup("warfarin", "aspirin")
        assert "major" in result["severity"].values

    def test_lookup_unknown_pair_returns_empty(self, db: DrugInteractionDB) -> None:
        """lookup() returns an empty DataFrame for an unknown drug pair."""
        result = db.lookup("penicillin", "acetaminophen")
        assert result.empty, "No interaction expected for unknown pair"

    def test_lookup_without_loaded_data_raises(self, db_empty: DrugInteractionDB) -> None:
        """lookup() raises RuntimeError when no data has been loaded."""
        with pytest.raises(RuntimeError, match="No data loaded"):
            db_empty.lookup("warfarin", "aspirin")


# ---------------------------------------------------------------------------
# Test 2: Severity filtering
# ---------------------------------------------------------------------------

class TestSeverityFilter:
    def test_filter_major_returns_only_major(self, db: DrugInteractionDB) -> None:
        """filter_by_severity('major') returns only rows with severity='major'."""
        result = db.filter_by_severity("major")
        assert not result.empty
        assert set(result["severity"].unique()) == {"major"}

    def test_filter_contraindicated(self, db: DrugInteractionDB) -> None:
        """filter_by_severity('contraindicated') returns contraindicated rows."""
        result = db.filter_by_severity("contraindicated")
        assert not result.empty
        assert all(result["severity"] == "contraindicated")

    def test_filter_minimum_threshold(self, db: DrugInteractionDB) -> None:
        """filter_by_severity with minimum=True returns major + contraindicated."""
        result = db.filter_by_severity("major", minimum=True)
        allowed = {"major", "contraindicated"}
        assert set(result["severity"].unique()).issubset(allowed)
        # Must include both severity levels present in sample
        assert "major" in result["severity"].values
        assert "contraindicated" in result["severity"].values

    def test_filter_minor_minimum_returns_all(self, db: DrugInteractionDB) -> None:
        """filter_by_severity('minor', minimum=True) returns all rows."""
        all_rows = db.filter_by_severity("minor", minimum=True)
        full_df = db._df
        assert len(all_rows) == len(full_df)

    def test_filter_unknown_severity_raises(self, db: DrugInteractionDB) -> None:
        """filter_by_severity raises ValueError for an unrecognised severity."""
        with pytest.raises(ValueError, match="Unknown severity"):
            db.filter_by_severity("extreme")

    def test_filter_case_insensitive(self, db: DrugInteractionDB) -> None:
        """filter_by_severity accepts severity in any case."""
        result_lower = db.filter_by_severity("major")
        result_upper = db.filter_by_severity("MAJOR")
        assert len(result_lower) == len(result_upper)


# ---------------------------------------------------------------------------
# Test 3: Bidirectional search
# ---------------------------------------------------------------------------

class TestBidirectionalSearch:
    def test_lookup_ab_equals_lookup_ba(self, db: DrugInteractionDB) -> None:
        """lookup(A, B) and lookup(B, A) return the same rows."""
        result_ab = db.lookup("warfarin", "aspirin")
        result_ba = db.lookup("aspirin", "warfarin")
        assert len(result_ab) == len(result_ba), (
            "Bidirectional lookup should return the same number of rows"
        )

    def test_bidirectional_simvastatin_amiodarone(self, db: DrugInteractionDB) -> None:
        """Reversed pair (amiodarone, simvastatin) finds the same interaction."""
        result = db.lookup("amiodarone", "simvastatin")
        assert not result.empty

    def test_search_drug_both_positions(self, db: DrugInteractionDB) -> None:
        """search_drug() finds a drug whether it appears in drug_a or drug_b."""
        # amiodarone appears as drug_b in simvastatin/amiodarone and digoxin/amiodarone
        result = db.search_drug("amiodarone")
        assert len(result) >= 2, (
            "amiodarone should appear in at least 2 interactions"
        )


# ---------------------------------------------------------------------------
# Test 4: Multi-drug interaction check
# ---------------------------------------------------------------------------

class TestMultiDrugCheck:
    def test_multi_drug_finds_interactions(self, db: DrugInteractionDB) -> None:
        """check_multi_drug returns interactions for a list with known pairs."""
        result = db.check_multi_drug(["warfarin", "aspirin", "clopidogrel"])
        assert not result.empty

    def test_multi_drug_less_than_two_raises(self, db: DrugInteractionDB) -> None:
        """check_multi_drug raises ValueError when fewer than 2 drugs given."""
        with pytest.raises(ValueError, match="at least two drugs"):
            db.check_multi_drug(["warfarin"])

    def test_multi_drug_no_interactions_returns_empty(self, db: DrugInteractionDB) -> None:
        """check_multi_drug returns empty DataFrame for drugs with no interactions."""
        result = db.check_multi_drug(["penicillin", "acetaminophen"])
        assert result.empty

    def test_multi_drug_queried_pair_column_present(self, db: DrugInteractionDB) -> None:
        """check_multi_drug result includes a queried_pair column."""
        result = db.check_multi_drug(["warfarin", "aspirin"])
        assert "queried_pair" in result.columns


# ---------------------------------------------------------------------------
# Test 5: Edge cases — unknown drug and self-interaction
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_unknown_drug_lookup_empty(self, db: DrugInteractionDB) -> None:
        """lookup() with an entirely unknown drug returns empty DataFrame."""
        result = db.lookup("unknowndrug123", "aspirin")
        assert result.empty

    def test_self_interaction_returns_empty(self, db: DrugInteractionDB) -> None:
        """lookup(drug, drug) returns empty since no self-interaction is recorded."""
        result = db.lookup("warfarin", "warfarin")
        assert result.empty

    def test_case_insensitive_lookup(self, db: DrugInteractionDB) -> None:
        """lookup() matches regardless of input capitalisation."""
        result_upper = db.lookup("WARFARIN", "ASPIRIN")
        result_mixed = db.lookup("Warfarin", "Aspirin")
        result_lower = db.lookup("warfarin", "aspirin")
        assert len(result_upper) == len(result_lower)
        assert len(result_mixed) == len(result_lower)

    def test_whitespace_in_drug_name_stripped(self, db: DrugInteractionDB) -> None:
        """lookup() strips surrounding whitespace from drug names."""
        result_padded = db.lookup("  warfarin  ", "  aspirin  ")
        result_clean = db.lookup("warfarin", "aspirin")
        assert len(result_padded) == len(result_clean)

    def test_validate_empty_dataframe_raises(self) -> None:
        """validate() raises ValueError on empty DataFrame."""
        instance = DrugInteractionDB()
        with pytest.raises(ValueError, match="empty"):
            instance.validate(pd.DataFrame())

    def test_validate_missing_required_columns_raises(self) -> None:
        """validate() raises ValueError when required columns are absent."""
        instance = DrugInteractionDB()
        bad_df = pd.DataFrame({"drug_a": ["aspirin"], "drug_b": ["warfarin"]})
        with pytest.raises(ValueError, match="Missing required columns"):
            instance.validate(bad_df)

    def test_preprocess_deduplicates(self) -> None:
        """preprocess() removes exact duplicate rows."""
        instance = DrugInteractionDB()
        base_row = {
            "drug_a": "Warfarin",
            "drug_b": "Aspirin",
            "interaction_type": "pharmacodynamic",
            "severity": "major",
        }
        df_with_dup = pd.DataFrame([base_row, base_row])
        result = instance.preprocess(df_with_dup)
        assert len(result) == 1, "Duplicate row should have been removed"

    def test_preprocess_does_not_mutate_input(self) -> None:
        """preprocess() returns a new DataFrame; the input is unchanged."""
        instance = DrugInteractionDB()
        original = pd.DataFrame([{
            "drug_a": "Warfarin",
            "drug_b": "Aspirin",
            "interaction_type": "PD",
            "severity": "Major",
        }])
        original_severity = original["severity"].iloc[0]
        instance.preprocess(original)
        assert original["severity"].iloc[0] == original_severity, (
            "preprocess() must not mutate the input DataFrame"
        )


# ---------------------------------------------------------------------------
# Test 6: list_drugs and get_contraindicated helpers
# ---------------------------------------------------------------------------

class TestHelpers:
    def test_list_drugs_returns_sorted_list(self, db: DrugInteractionDB) -> None:
        """list_drugs() returns a sorted list of unique drug names."""
        drugs = db.list_drugs()
        assert isinstance(drugs, list)
        assert drugs == sorted(drugs)

    def test_list_drugs_all_lowercase(self, db: DrugInteractionDB) -> None:
        """list_drugs() returns lowercase drug names."""
        drugs = db.list_drugs()
        assert all(d == d.lower() for d in drugs)

    def test_get_contraindicated_subset(self, db: DrugInteractionDB) -> None:
        """get_contraindicated() returns only contraindicated rows."""
        result = db.get_contraindicated()
        assert not result.empty
        assert all(result["severity"] == "contraindicated")

    def test_analyze_returns_required_keys(self, db: DrugInteractionDB) -> None:
        """analyze() result includes total_records, columns, and severity_counts."""
        result = db.analyze(db._df)
        assert "total_records" in result
        assert "columns" in result
        assert "severity_counts" in result

    def test_to_dataframe_has_metric_value_columns(self, db: DrugInteractionDB) -> None:
        """to_dataframe() produces a DataFrame with 'metric' and 'value' columns."""
        analysis = db.analyze(db._df)
        flat = db.to_dataframe(analysis)
        assert "metric" in flat.columns
        assert "value" in flat.columns
