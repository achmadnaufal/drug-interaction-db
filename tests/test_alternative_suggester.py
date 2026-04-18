"""
Tests for src/alternative_suggester.py.

Run with:
    pytest tests/test_alternative_suggester.py -v
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.alternative_suggester import (
    THERAPEUTIC_CLASSES,
    AlternativeSuggestion,
    list_classes,
    suggest_alternatives,
    suggestions_to_dataframe,
)
from src.main import DrugInteractionDB


SAMPLE_CSV = str(Path(__file__).parent.parent / "demo" / "sample_data.csv")


@pytest.fixture(scope="module")
def loaded_db() -> DrugInteractionDB:
    """Return a DrugInteractionDB loaded with the sample dataset."""
    db = DrugInteractionDB()
    db.load_data(SAMPLE_CSV)
    return db


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestHappyPath:
    def test_returns_tuple(self, loaded_db: DrugInteractionDB) -> None:
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "fluconazole"],
            db=loaded_db,
        )
        assert isinstance(result, tuple)

    def test_each_item_is_suggestion(self, loaded_db: DrugInteractionDB) -> None:
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "fluconazole"],
            db=loaded_db,
        )
        for item in result:
            assert isinstance(item, AlternativeSuggestion)

    def test_results_sorted_ascending_by_score(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "fluconazole", "tacrolimus"],
            db=loaded_db,
            max_results=10,
        )
        scores = [s.total_score for s in result]
        assert scores == sorted(scores)

    def test_excludes_target_from_candidates(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "fluconazole"],
            db=loaded_db,
        )
        names = [s.candidate for s in result]
        assert "fluconazole" not in names

    def test_case_insensitive_target(self, loaded_db: DrugInteractionDB) -> None:
        upper = suggest_alternatives(
            target="FLUCONAZOLE",
            regimen=["warfarin"],
            db=loaded_db,
        )
        lower = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin"],
            db=loaded_db,
        )
        assert [s.candidate for s in upper] == [s.candidate for s in lower]

    def test_case_insensitive_regimen(self, loaded_db: DrugInteractionDB) -> None:
        mixed = suggest_alternatives(
            target="fluconazole",
            regimen=["Warfarin", "FLUCONAZOLE"],
            db=loaded_db,
        )
        assert isinstance(mixed, tuple)
        # Target should not leak into candidates regardless of case.
        assert all(s.candidate != "fluconazole" for s in mixed)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_unknown_drug_class_returns_empty(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        result = suggest_alternatives(
            target="zzz_not_a_drug",
            regimen=["warfarin"],
            db=loaded_db,
        )
        assert result == ()

    def test_single_regimen_member_still_scores(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        # Only target itself in regimen — no "other" drugs to clash with.
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["fluconazole"],
            db=loaded_db,
        )
        # All candidates should score 0 because nothing else in regimen.
        assert all(s.total_score == 0.0 for s in result)
        assert all(s.highest_severity == "none" for s in result)

    def test_duplicate_regimen_entries_deduplicated(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        single = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "fluconazole"],
            db=loaded_db,
        )
        dup = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "warfarin", "fluconazole", "fluconazole"],
            db=loaded_db,
        )
        # Same candidates, same scores.
        assert [(s.candidate, s.total_score) for s in single] == [
            (s.candidate, s.total_score) for s in dup
        ]

    def test_max_results_caps_output(self, loaded_db: DrugInteractionDB) -> None:
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin"],
            db=loaded_db,
            max_results=2,
        )
        assert len(result) <= 2

    def test_is_safer_flag_reflects_baseline(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        # Fluconazole + warfarin has a known major interaction, so
        # fluconazole's baseline score is > 0.  Candidates that don't
        # appear in the DB at all should score 0 and therefore be safer.
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin", "fluconazole"],
            db=loaded_db,
        )
        zero_scored = [s for s in result if s.total_score == 0.0]
        assert all(s.is_safer for s in zero_scored)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:
    def test_empty_target_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError):
            suggest_alternatives(target="", regimen=["warfarin"], db=loaded_db)

    def test_whitespace_target_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError):
            suggest_alternatives(target="   ", regimen=["warfarin"], db=loaded_db)

    def test_non_string_target_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(TypeError):
            suggest_alternatives(target=123, regimen=["warfarin"], db=loaded_db)  # type: ignore[arg-type]

    def test_empty_regimen_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(ValueError):
            suggest_alternatives(target="fluconazole", regimen=[], db=loaded_db)

    def test_non_list_regimen_raises(self, loaded_db: DrugInteractionDB) -> None:
        with pytest.raises(TypeError):
            suggest_alternatives(
                target="fluconazole",
                regimen="warfarin",  # type: ignore[arg-type]
                db=loaded_db,
            )

    def test_blank_regimen_entry_raises(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        with pytest.raises(ValueError):
            suggest_alternatives(
                target="fluconazole",
                regimen=["warfarin", ""],
                db=loaded_db,
            )

    def test_non_db_raises(self) -> None:
        with pytest.raises(TypeError):
            suggest_alternatives(
                target="fluconazole",
                regimen=["warfarin"],
                db="not a db",  # type: ignore[arg-type]
            )

    def test_zero_max_results_raises(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        with pytest.raises(ValueError):
            suggest_alternatives(
                target="fluconazole",
                regimen=["warfarin"],
                db=loaded_db,
                max_results=0,
            )

    def test_bool_max_results_raises(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        # bool is a subclass of int but should be rejected explicitly.
        with pytest.raises(TypeError):
            suggest_alternatives(
                target="fluconazole",
                regimen=["warfarin"],
                db=loaded_db,
                max_results=True,  # type: ignore[arg-type]
            )


# ---------------------------------------------------------------------------
# DataFrame flattening
# ---------------------------------------------------------------------------

class TestSuggestionsToDataFrame:
    def test_empty_input_returns_empty_df(self) -> None:
        df = suggestions_to_dataframe(())
        assert isinstance(df, pd.DataFrame)
        assert df.empty
        assert list(df.columns) == [
            "candidate",
            "therapeutic_class",
            "total_score",
            "highest_severity",
            "interacting_drugs",
            "is_safer",
        ]

    def test_populated_input_returns_rows(
        self, loaded_db: DrugInteractionDB
    ) -> None:
        result = suggest_alternatives(
            target="fluconazole",
            regimen=["warfarin"],
            db=loaded_db,
        )
        df = suggestions_to_dataframe(result)
        assert len(df) == len(result)
        assert list(df.columns) == [
            "candidate",
            "therapeutic_class",
            "total_score",
            "highest_severity",
            "interacting_drugs",
            "is_safer",
        ]

    def test_non_tuple_input_raises(self) -> None:
        with pytest.raises(TypeError):
            suggestions_to_dataframe([])  # type: ignore[arg-type]

    def test_bad_item_in_tuple_raises(self) -> None:
        with pytest.raises(TypeError):
            suggestions_to_dataframe(("not a suggestion",))  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Class registry
# ---------------------------------------------------------------------------

class TestClassRegistry:
    def test_list_classes_returns_tuple(self) -> None:
        classes = list_classes()
        assert isinstance(classes, tuple)
        assert len(classes) > 0

    def test_all_classes_have_members(self) -> None:
        for name, members in THERAPEUTIC_CLASSES.items():
            assert isinstance(name, str) and name == name.lower()
            assert isinstance(members, tuple)
            assert len(members) >= 2

    def test_all_members_are_lowercase(self) -> None:
        for members in THERAPEUTIC_CLASSES.values():
            for m in members:
                assert m == m.lower()
