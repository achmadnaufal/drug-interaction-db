"""
Tests for src/polypharmacy_risk_scorer.py.

Run with:
    pytest tests/test_polypharmacy_risk_scorer.py -v
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.main import DrugInteractionDB
from src.polypharmacy_risk_scorer import (
    HIGH_BURDEN_THRESHOLD,
    SEVERITY_WEIGHTS,
    PairScore,
    RiskReport,
    report_to_dataframe,
    score_regimen,
)

SAMPLE_CSV = str(Path(__file__).parent.parent / "demo" / "sample_data.csv")


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def loaded_db() -> DrugInteractionDB:
    """Return a DrugInteractionDB loaded with the sample dataset."""
    db = DrugInteractionDB()
    db.load_data(SAMPLE_CSV)
    return db


# ---------------------------------------------------------------------------
# 1. Happy path — known interacting pair
# ---------------------------------------------------------------------------

class TestHappyPath:
    def test_known_pair_returns_nonzero_score(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "aspirin"], loaded_db)
        assert report.total_score > 0

    def test_known_pair_interacting_pairs_count(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "aspirin"], loaded_db)
        assert report.interacting_pairs == 1

    def test_known_pair_severity_is_major(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "aspirin"], loaded_db)
        assert report.pair_scores[0].severity == "major"

    def test_known_pair_score_matches_weight(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "aspirin"], loaded_db)
        ps = report.pair_scores[0]
        assert ps.pair_score == pytest.approx(SEVERITY_WEIGHTS["major"] * ps.interaction_count)

    def test_three_drug_regimen_aggregates_correctly(self, loaded_db: DrugInteractionDB) -> None:
        # warfarin-aspirin (major) and warfarin-fluconazole (major) should both fire
        report = score_regimen(["warfarin", "aspirin", "fluconazole"], loaded_db)
        assert report.interacting_pairs >= 2
        assert report.total_score >= SEVERITY_WEIGHTS["major"] * 2


# ---------------------------------------------------------------------------
# 2. Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_single_drug_returns_zero_score(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin"], loaded_db)
        assert report.total_score == 0.0
        assert report.interacting_pairs == 0
        assert report.total_pairs == 0
        assert report.risk_level == "low"

    def test_two_unknown_drugs_no_interactions(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["unknowndrug_x", "unknowndrug_y"], loaded_db)
        assert report.total_score == 0.0
        assert report.interacting_pairs == 0
        assert report.risk_level == "low"

    def test_duplicate_drugs_deduplicated(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "warfarin", "aspirin"], loaded_db)
        # After dedup only warfarin + aspirin remain — same as two-drug call
        report2 = score_regimen(["warfarin", "aspirin"], loaded_db)
        assert report.total_score == report2.total_score

    def test_contraindicated_pair_flags_has_contraindicated(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["ssris", "maois"], loaded_db)
        assert report.has_contraindicated is True

    def test_contraindicated_pair_risk_level_is_critical(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["ssris", "maois"], loaded_db)
        assert report.risk_level == "critical"


# ---------------------------------------------------------------------------
# 3. Polypharmacy threshold
# ---------------------------------------------------------------------------

class TestPolypharmacyFlag:
    def test_five_drugs_flags_polypharmacy(self, loaded_db: DrugInteractionDB) -> None:
        drugs = ["warfarin", "aspirin", "fluconazole", "simvastatin", "amiodarone"]
        report = score_regimen(drugs, loaded_db)
        assert report.polypharmacy_flagged is True

    def test_four_drugs_does_not_flag_polypharmacy(self, loaded_db: DrugInteractionDB) -> None:
        drugs = ["warfarin", "aspirin", "fluconazole", "simvastatin"]
        report = score_regimen(drugs, loaded_db)
        assert report.polypharmacy_flagged is False

    def test_custom_polypharmacy_count(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "aspirin"], loaded_db, polypharmacy_count=2)
        assert report.polypharmacy_flagged is True


# ---------------------------------------------------------------------------
# 4. Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:
    def test_empty_list_raises_value_error(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            score_regimen([], loaded_db)

    def test_non_list_raises_type_error(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(TypeError, match="must be a list"):
            score_regimen("warfarin", loaded_db)  # type: ignore[arg-type]

    def test_non_db_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="DrugInteractionDB"):
            score_regimen(["warfarin"], "not_a_db")  # type: ignore[arg-type]

    def test_blank_drug_name_raises_value_error(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError):
            score_regimen(["warfarin", "   "], loaded_db)

    def test_non_string_drug_name_raises_type_error(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(TypeError):
            score_regimen([123, "aspirin"], loaded_db)  # type: ignore[list-item]


# ---------------------------------------------------------------------------
# 5. Determinism
# ---------------------------------------------------------------------------

class TestDeterminism:
    def test_same_input_same_score(self, loaded_db: DrugInteractionDB) -> None:
        drugs = ["warfarin", "aspirin", "fluconazole"]
        r1 = score_regimen(drugs, loaded_db)
        r2 = score_regimen(drugs, loaded_db)
        assert r1.total_score == r2.total_score
        assert r1.risk_level == r2.risk_level

    def test_order_independent_score(self, loaded_db: DrugInteractionDB) -> None:
        drugs_fwd = ["warfarin", "aspirin"]
        drugs_rev = ["aspirin", "warfarin"]
        r1 = score_regimen(drugs_fwd, loaded_db)
        r2 = score_regimen(drugs_rev, loaded_db)
        assert r1.total_score == r2.total_score


# ---------------------------------------------------------------------------
# 6. report_to_dataframe helper
# ---------------------------------------------------------------------------

class TestReportToDataframe:
    def test_empty_report_returns_empty_df(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["unknowndrug_x", "unknowndrug_y"], loaded_db)
        df = report_to_dataframe(report)
        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_non_empty_report_has_correct_columns(self, loaded_db: DrugInteractionDB) -> None:
        report = score_regimen(["warfarin", "aspirin"], loaded_db)
        df = report_to_dataframe(report)
        expected_cols = {"drug_a", "drug_b", "severity", "weight", "interaction_count", "pair_score"}
        assert expected_cols.issubset(set(df.columns))

    def test_invalid_type_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="RiskReport"):
            report_to_dataframe({"not": "a report"})  # type: ignore[arg-type]

    def test_df_sorted_by_pair_score_descending(self, loaded_db: DrugInteractionDB) -> None:
        drugs = ["warfarin", "aspirin", "fluconazole", "simvastatin", "amiodarone"]
        report = score_regimen(drugs, loaded_db)
        df = report_to_dataframe(report)
        if len(df) > 1:
            assert df["pair_score"].iloc[0] >= df["pair_score"].iloc[-1]


# ---------------------------------------------------------------------------
# 7. Parametrized severity weights
# ---------------------------------------------------------------------------

class TestCustomWeights:
    @pytest.mark.parametrize("severity,expected_weight", [
        ("minor", 1.0),
        ("moderate", 3.0),
        ("major", 7.0),
        ("contraindicated", 15.0),
    ])
    def test_default_weight_values(self, severity: str, expected_weight: float) -> None:
        assert SEVERITY_WEIGHTS[severity] == pytest.approx(expected_weight)

    def test_custom_weights_applied(self, loaded_db: DrugInteractionDB) -> None:
        custom = {"minor": 1.0, "moderate": 1.0, "major": 1.0, "contraindicated": 1.0}
        report = score_regimen(["warfarin", "aspirin"], loaded_db, severity_weights=custom)
        # All weights are 1.0, so pair_score == interaction_count
        ps = report.pair_scores[0]
        assert ps.pair_score == pytest.approx(float(ps.interaction_count))
