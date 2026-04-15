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
  - Type errors on invalid inputs
  - has_interaction helper
  - get_high_risk_drugs helper
  - load_data file handling
"""

import pytest
import pandas as pd
import sys
from pathlib import Path

# Ensure project root is importable regardless of where pytest is invoked
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.main import DrugInteractionDB, SEVERITY_ORDER, REQUIRED_COLUMNS


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
# Module-level constants
# ---------------------------------------------------------------------------

class TestConstants:
    def test_severity_order_has_four_levels(self) -> None:
        """SEVERITY_ORDER contains exactly four recognised severity levels."""
        assert len(SEVERITY_ORDER) == 4

    def test_severity_order_sequence(self) -> None:
        """SEVERITY_ORDER increases from minor to contraindicated."""
        assert SEVERITY_ORDER == ["minor", "moderate", "major", "contraindicated"]

    def test_required_columns_is_frozenset(self) -> None:
        """REQUIRED_COLUMNS is a frozenset (immutable)."""
        assert isinstance(REQUIRED_COLUMNS, frozenset)

    def test_required_columns_contains_minimum_fields(self) -> None:
        """REQUIRED_COLUMNS includes the four mandatory fields."""
        assert {"drug_a", "drug_b", "interaction_type", "severity"}.issubset(
            REQUIRED_COLUMNS
        )


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

    def test_lookup_returns_copy_not_reference(self, db: DrugInteractionDB) -> None:
        """lookup() returns a copy so callers cannot mutate internal state."""
        result = db.lookup("warfarin", "aspirin")
        original_len = len(db._df)
        result.drop(result.index, inplace=True)
        assert len(db._df) == original_len, "Internal DataFrame should not be mutated"


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

    def test_filter_minor_exact(self, db: DrugInteractionDB) -> None:
        """filter_by_severity('minor') returns only minor-severity rows."""
        result = db.filter_by_severity("minor")
        assert not result.empty
        assert all(result["severity"] == "minor")

    def test_filter_moderate_exact(self, db: DrugInteractionDB) -> None:
        """filter_by_severity('moderate') returns only moderate rows."""
        result = db.filter_by_severity("moderate")
        assert not result.empty
        assert all(result["severity"] == "moderate")

    def test_filter_returns_copy(self, db: DrugInteractionDB) -> None:
        """filter_by_severity() returns a copy; mutating it does not affect internal state."""
        result = db.filter_by_severity("major")
        original_len = len(db._df)
        result.drop(result.index, inplace=True)
        assert len(db._df) == original_len

    def test_filter_severity_type_error(self, db: DrugInteractionDB) -> None:
        """filter_by_severity raises TypeError when severity is not a string."""
        with pytest.raises(TypeError):
            db.filter_by_severity(3)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Test 3: Bidirectional search
# ---------------------------------------------------------------------------

class TestBidirectionalSearch:
    def test_lookup_ab_equals_lookup_ba(self, db: DrugInteractionDB) -> None:
        """lookup(A, B) and lookup(B, A) return the same number of rows."""
        result_ab = db.lookup("warfarin", "aspirin")
        result_ba = db.lookup("aspirin", "warfarin")
        assert len(result_ab) == len(result_ba)

    def test_bidirectional_simvastatin_amiodarone(self, db: DrugInteractionDB) -> None:
        """Reversed pair (amiodarone, simvastatin) finds the same interaction."""
        result = db.lookup("amiodarone", "simvastatin")
        assert not result.empty

    def test_search_drug_both_positions(self, db: DrugInteractionDB) -> None:
        """search_drug() finds a drug whether it appears in drug_a or drug_b."""
        # amiodarone appears as drug_b in simvastatin/amiodarone and digoxin/amiodarone
        result = db.search_drug("amiodarone")
        assert len(result) >= 2

    def test_search_drug_case_insensitive(self, db: DrugInteractionDB) -> None:
        """search_drug() matches regardless of input capitalisation."""
        result_lower = db.search_drug("amiodarone")
        result_upper = db.search_drug("AMIODARONE")
        assert len(result_lower) == len(result_upper)

    def test_search_drug_empty_raises(self, db: DrugInteractionDB) -> None:
        """search_drug() raises ValueError for an empty drug name."""
        with pytest.raises(ValueError, match="must not be empty"):
            db.search_drug("")

    def test_search_drug_type_error(self, db: DrugInteractionDB) -> None:
        """search_drug() raises TypeError when drug is not a string."""
        with pytest.raises(TypeError):
            db.search_drug(None)  # type: ignore[arg-type]

    def test_search_unknown_drug_returns_empty(self, db: DrugInteractionDB) -> None:
        """search_drug() returns empty DataFrame for a drug not in the dataset."""
        result = db.search_drug("unknownxyz123")
        assert result.empty


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

    def test_multi_drug_empty_list_raises(self, db: DrugInteractionDB) -> None:
        """check_multi_drug raises ValueError for an empty list."""
        with pytest.raises(ValueError, match="at least two drugs"):
            db.check_multi_drug([])

    def test_multi_drug_non_list_raises(self, db: DrugInteractionDB) -> None:
        """check_multi_drug raises TypeError when drugs is not a list."""
        with pytest.raises(TypeError):
            db.check_multi_drug("warfarin,aspirin")  # type: ignore[arg-type]

    def test_multi_drug_no_interactions_returns_empty(self, db: DrugInteractionDB) -> None:
        """check_multi_drug returns empty DataFrame for drugs with no interactions."""
        result = db.check_multi_drug(["penicillin", "acetaminophen"])
        assert result.empty

    def test_multi_drug_queried_pair_column_present(self, db: DrugInteractionDB) -> None:
        """check_multi_drug result includes a queried_pair column."""
        result = db.check_multi_drug(["warfarin", "aspirin"])
        assert "queried_pair" in result.columns

    def test_multi_drug_three_way_all_pairs_checked(self, db: DrugInteractionDB) -> None:
        """check_multi_drug examines all n*(n-1)/2 pairs for a 3-drug list."""
        # simvastatin+amiodarone and digoxin+amiodarone are both in fixture
        result = db.check_multi_drug(["simvastatin", "amiodarone", "digoxin"])
        pairs = set(result["queried_pair"].unique())
        assert len(pairs) >= 2

    def test_multi_drug_empty_drug_name_raises(self, db: DrugInteractionDB) -> None:
        """check_multi_drug raises ValueError if any entry is an empty string."""
        with pytest.raises(ValueError):
            db.check_multi_drug(["warfarin", ""])


# ---------------------------------------------------------------------------
# Test 5: Edge cases — unknown drug, self-interaction, and type safety
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

    def test_lookup_empty_drug_a_raises(self, db: DrugInteractionDB) -> None:
        """lookup() raises ValueError when drug_a is empty."""
        with pytest.raises(ValueError, match="must not be empty"):
            db.lookup("", "aspirin")

    def test_lookup_empty_drug_b_raises(self, db: DrugInteractionDB) -> None:
        """lookup() raises ValueError when drug_b is empty."""
        with pytest.raises(ValueError, match="must not be empty"):
            db.lookup("warfarin", "")

    def test_lookup_whitespace_only_drug_raises(self, db: DrugInteractionDB) -> None:
        """lookup() raises ValueError when drug name is whitespace-only."""
        with pytest.raises(ValueError, match="must not be empty"):
            db.lookup("   ", "aspirin")

    def test_lookup_none_drug_raises(self, db: DrugInteractionDB) -> None:
        """lookup() raises TypeError when drug name is None."""
        with pytest.raises(TypeError):
            db.lookup(None, "aspirin")  # type: ignore[arg-type]

    def test_lookup_int_drug_raises(self, db: DrugInteractionDB) -> None:
        """lookup() raises TypeError when drug name is an integer."""
        with pytest.raises(TypeError):
            db.lookup(42, "aspirin")  # type: ignore[arg-type]

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

    def test_validate_non_dataframe_raises(self) -> None:
        """validate() raises TypeError when given a non-DataFrame argument."""
        instance = DrugInteractionDB()
        with pytest.raises(TypeError):
            instance.validate({"drug_a": ["aspirin"]})  # type: ignore[arg-type]

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
        assert original["severity"].iloc[0] == original_severity

    def test_preprocess_normalises_case(self) -> None:
        """preprocess() lowercases drug names and severity."""
        instance = DrugInteractionDB()
        df = pd.DataFrame([{
            "drug_a": "WARFARIN",
            "drug_b": "ASPIRIN",
            "interaction_type": "PHARMACODYNAMIC",
            "severity": "MAJOR",
        }])
        result = instance.preprocess(df)
        assert result["drug_a"].iloc[0] == "warfarin"
        assert result["drug_b"].iloc[0] == "aspirin"
        assert result["severity"].iloc[0] == "major"

    def test_preprocess_strips_whitespace(self) -> None:
        """preprocess() strips leading/trailing whitespace from string values."""
        instance = DrugInteractionDB()
        df = pd.DataFrame([{
            "drug_a": "  warfarin  ",
            "drug_b": "  aspirin  ",
            "interaction_type": " pharmacodynamic ",
            "severity": " major ",
        }])
        result = instance.preprocess(df)
        assert result["drug_a"].iloc[0] == "warfarin"
        assert result["severity"].iloc[0] == "major"

    def test_preprocess_non_dataframe_raises(self) -> None:
        """preprocess() raises TypeError when given a non-DataFrame argument."""
        instance = DrugInteractionDB()
        with pytest.raises(TypeError):
            instance.preprocess([{"drug_a": "a"}])  # type: ignore[arg-type]


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

    def test_list_drugs_no_duplicates(self, db: DrugInteractionDB) -> None:
        """list_drugs() contains no duplicate entries."""
        drugs = db.list_drugs()
        assert len(drugs) == len(set(drugs))

    def test_list_drugs_without_loaded_data_raises(self, db_empty: DrugInteractionDB) -> None:
        """list_drugs() raises RuntimeError when no data is loaded."""
        with pytest.raises(RuntimeError, match="No data loaded"):
            db_empty.list_drugs()

    def test_get_contraindicated_subset(self, db: DrugInteractionDB) -> None:
        """get_contraindicated() returns only contraindicated rows."""
        result = db.get_contraindicated()
        assert not result.empty
        assert all(result["severity"] == "contraindicated")

    def test_get_contraindicated_without_data_raises(self, db_empty: DrugInteractionDB) -> None:
        """get_contraindicated() raises RuntimeError when no data is loaded."""
        with pytest.raises(RuntimeError, match="No data loaded"):
            db_empty.get_contraindicated()

    def test_analyze_returns_required_keys(self, db: DrugInteractionDB) -> None:
        """analyze() result includes total_records, columns, and severity_counts."""
        result = db.analyze(db._df)
        assert "total_records" in result
        assert "columns" in result
        assert "severity_counts" in result

    def test_analyze_total_records_correct(self, db: DrugInteractionDB) -> None:
        """analyze() total_records matches the number of rows in the DataFrame."""
        result = db.analyze(db._df)
        assert result["total_records"] == len(db._df)

    def test_analyze_severity_counts_sum(self, db: DrugInteractionDB) -> None:
        """Sum of severity_counts equals total_records."""
        result = db.analyze(db._df)
        counts_sum = sum(result["severity_counts"].values())
        assert counts_sum == result["total_records"]

    def test_to_dataframe_has_metric_value_columns(self, db: DrugInteractionDB) -> None:
        """to_dataframe() produces a DataFrame with 'metric' and 'value' columns."""
        analysis = db.analyze(db._df)
        flat = db.to_dataframe(analysis)
        assert "metric" in flat.columns
        assert "value" in flat.columns

    def test_to_dataframe_non_dict_raises(self, db: DrugInteractionDB) -> None:
        """to_dataframe() raises TypeError for non-dict input."""
        with pytest.raises(TypeError):
            db.to_dataframe("not a dict")  # type: ignore[arg-type]

    def test_to_dataframe_nested_keys_dotted(self, db: DrugInteractionDB) -> None:
        """to_dataframe() expands nested dicts with dotted metric names."""
        analysis = db.analyze(db._df)
        flat = db.to_dataframe(analysis)
        dotted = [m for m in flat["metric"] if "." in m]
        assert len(dotted) > 0, "Expected at least one dotted metric from severity_counts"


# ---------------------------------------------------------------------------
# Test 7: has_interaction helper
# ---------------------------------------------------------------------------

class TestHasInteraction:
    def test_has_interaction_true_for_known_pair(self, db: DrugInteractionDB) -> None:
        """has_interaction() returns True for a known drug pair."""
        assert db.has_interaction("warfarin", "aspirin") is True

    def test_has_interaction_false_for_unknown_pair(self, db: DrugInteractionDB) -> None:
        """has_interaction() returns False for drugs with no recorded interaction."""
        assert db.has_interaction("penicillin", "acetaminophen") is False

    def test_has_interaction_bidirectional(self, db: DrugInteractionDB) -> None:
        """has_interaction() is True regardless of argument order."""
        assert db.has_interaction("aspirin", "warfarin") is True

    def test_has_interaction_empty_name_raises(self, db: DrugInteractionDB) -> None:
        """has_interaction() raises ValueError for an empty drug name."""
        with pytest.raises(ValueError):
            db.has_interaction("", "aspirin")


# ---------------------------------------------------------------------------
# Test 8: get_high_risk_drugs helper
# ---------------------------------------------------------------------------

class TestGetHighRiskDrugs:
    def test_returns_list(self, db: DrugInteractionDB) -> None:
        """get_high_risk_drugs() returns a list."""
        result = db.get_high_risk_drugs()
        assert isinstance(result, list)

    def test_result_is_sorted(self, db: DrugInteractionDB) -> None:
        """get_high_risk_drugs() returns drugs in sorted order."""
        result = db.get_high_risk_drugs()
        assert result == sorted(result)

    def test_result_all_lowercase(self, db: DrugInteractionDB) -> None:
        """get_high_risk_drugs() returns lowercase drug names."""
        result = db.get_high_risk_drugs()
        assert all(d == d.lower() for d in result)

    def test_warfarin_in_high_risk(self, db: DrugInteractionDB) -> None:
        """warfarin appears in high-risk drugs since it has a 'major' interaction."""
        result = db.get_high_risk_drugs()
        assert "warfarin" in result

    def test_no_minor_only_drug_included(self, db: DrugInteractionDB) -> None:
        """Drugs that only appear in 'minor' interactions are not in high-risk list."""
        # phenytoin/folic_acid is the only 'minor' pair in SAMPLE_ROWS
        result = db.get_high_risk_drugs()
        # folic_acid only appears in the minor pair, so it must not be in high_risk
        assert "folic_acid" not in result

    def test_without_data_raises(self, db_empty: DrugInteractionDB) -> None:
        """get_high_risk_drugs() raises RuntimeError when no data is loaded."""
        with pytest.raises(RuntimeError, match="No data loaded"):
            db_empty.get_high_risk_drugs()


# ---------------------------------------------------------------------------
# Test 9: load_data error handling
# ---------------------------------------------------------------------------

class TestLoadData:
    def test_load_nonexistent_file_raises(self, db_empty: DrugInteractionDB) -> None:
        """load_data() raises FileNotFoundError for a missing file."""
        with pytest.raises(FileNotFoundError):
            db_empty.load_data("/tmp/does_not_exist_xyz.csv")

    def test_load_non_string_filepath_raises(self, db_empty: DrugInteractionDB) -> None:
        """load_data() raises TypeError when filepath is not a string."""
        with pytest.raises(TypeError):
            db_empty.load_data(123)  # type: ignore[arg-type]

    def test_load_sample_csv(self, tmp_path: Path) -> None:
        """load_data() successfully loads a valid CSV file."""
        csv_path = tmp_path / "test.csv"
        df = pd.DataFrame(SAMPLE_ROWS)
        df.to_csv(csv_path, index=False)

        instance = DrugInteractionDB()
        result = instance.load_data(str(csv_path))
        assert not result.empty
        assert "drug_a" in result.columns

    def test_load_sets_internal_df(self, tmp_path: Path) -> None:
        """load_data() populates _df so subsequent query methods work."""
        csv_path = tmp_path / "test.csv"
        df = pd.DataFrame(SAMPLE_ROWS)
        df.to_csv(csv_path, index=False)

        instance = DrugInteractionDB()
        instance.load_data(str(csv_path))
        assert instance._df is not None

    def test_load_empty_csv_raises(self, tmp_path: Path) -> None:
        """load_data() raises ValueError when the file has no data rows."""
        csv_path = tmp_path / "empty.csv"
        pd.DataFrame(columns=["drug_a", "drug_b", "interaction_type", "severity"]).to_csv(
            csv_path, index=False
        )
        instance = DrugInteractionDB()
        with pytest.raises(ValueError):
            instance.load_data(str(csv_path))

    def test_load_missing_column_raises(self, tmp_path: Path) -> None:
        """load_data() raises ValueError when required columns are missing."""
        csv_path = tmp_path / "bad.csv"
        pd.DataFrame({"drug_a": ["warfarin"], "drug_b": ["aspirin"]}).to_csv(
            csv_path, index=False
        )
        instance = DrugInteractionDB()
        with pytest.raises(ValueError, match="Missing required columns"):
            instance.load_data(str(csv_path))


# ---------------------------------------------------------------------------
# Test 10: DrugInteractionDB config and init
# ---------------------------------------------------------------------------

class TestInit:
    def test_default_config_is_empty_dict(self) -> None:
        """DrugInteractionDB() with no args has an empty config dict."""
        db = DrugInteractionDB()
        assert db.config == {}

    def test_custom_config_stored(self) -> None:
        """DrugInteractionDB(config=...) stores the provided config."""
        cfg = {"strict_validation": True}
        db = DrugInteractionDB(config=cfg)
        assert db.config == cfg

    def test_config_does_not_mutate_default(self) -> None:
        """Passing a config dict does not share state between instances."""
        db1 = DrugInteractionDB(config={"key": "a"})
        db2 = DrugInteractionDB(config={"key": "b"})
        assert db1.config["key"] != db2.config["key"]

    def test_initial_df_is_none(self) -> None:
        """_df is None before any data is loaded."""
        db = DrugInteractionDB()
        assert db._df is None


# ---------------------------------------------------------------------------
# Test 11: data_generator module
# ---------------------------------------------------------------------------

class TestDataGenerator:
    def test_generate_sample_returns_dataframe(self) -> None:
        """generate_sample() returns a pandas DataFrame."""
        from src.data_generator import generate_sample
        df = generate_sample()
        assert isinstance(df, pd.DataFrame)

    def test_generate_sample_default_20_rows(self) -> None:
        """generate_sample() returns 20 rows by default."""
        from src.data_generator import generate_sample
        df = generate_sample()
        assert len(df) == 20

    def test_generate_sample_respects_n(self) -> None:
        """generate_sample(n=5) returns at most 5 rows."""
        from src.data_generator import generate_sample
        df = generate_sample(n=5)
        assert len(df) == 5

    def test_generate_sample_has_required_columns(self) -> None:
        """generate_sample() output contains all columns required by DrugInteractionDB."""
        from src.data_generator import generate_sample
        df = generate_sample()
        for col in ("drug_a", "drug_b", "interaction_type", "severity"):
            assert col in df.columns, f"Missing required column: {col}"

    def test_generate_sample_severity_values_valid(self) -> None:
        """All severity values in generated data are recognised severity levels."""
        from src.data_generator import generate_sample, COLUMNS
        df = generate_sample()
        valid = {"minor", "moderate", "major", "contraindicated"}
        actual = set(df["severity"].str.lower().unique())
        assert actual.issubset(valid), f"Unexpected severity values: {actual - valid}"

    def test_generate_sample_no_empty_drug_names(self) -> None:
        """Generated drug names are non-empty strings."""
        from src.data_generator import generate_sample
        df = generate_sample()
        assert df["drug_a"].notna().all()
        assert df["drug_b"].notna().all()
        assert (df["drug_a"].str.strip() != "").all()
        assert (df["drug_b"].str.strip() != "").all()

    def test_generate_sample_loadable_by_db(self, tmp_path: Path) -> None:
        """Data from generate_sample() can be loaded and queried by DrugInteractionDB."""
        from src.data_generator import generate_sample
        df = generate_sample()
        csv_path = tmp_path / "generated.csv"
        df.to_csv(csv_path, index=False)

        instance = DrugInteractionDB()
        loaded = instance.load_data(str(csv_path))
        assert not loaded.empty
        assert len(instance.list_drugs()) > 0
