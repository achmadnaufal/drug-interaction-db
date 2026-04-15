"""
Drug-drug interaction database query and risk assessment tools.

This module provides the DrugInteractionDB class for loading, validating,
querying, and analyzing drug-drug interaction data from CSV or Excel files.

Author: github.com/achmadnaufal
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List


SEVERITY_ORDER: List[str] = ["minor", "moderate", "major", "contraindicated"]

REQUIRED_COLUMNS: frozenset = frozenset({"drug_a", "drug_b", "interaction_type", "severity"})


def _validate_drug_name(name: str, param: str = "drug") -> str:
    """
    Validate and normalise a drug name string.

    Args:
        name: Raw drug name provided by the caller.
        param: Parameter label used in error messages (e.g. ``"drug_a"``).

    Returns:
        Stripped, lowercased drug name ready for database matching.

    Raises:
        TypeError: If *name* is not a string.
        ValueError: If *name* is empty or whitespace-only after stripping.
    """
    if not isinstance(name, str):
        raise TypeError(
            f"'{param}' must be a string, got {type(name).__name__!r}."
        )
    cleaned = name.strip()
    if not cleaned:
        raise ValueError(
            f"'{param}' must not be empty or whitespace-only."
        )
    return cleaned.lower()


class DrugInteractionDB:
    """
    Drug interaction database query and risk assessment tool.

    Loads drug-drug interaction data from CSV or Excel files and provides
    methods for lookup, severity filtering, bidirectional search, and
    multi-drug interaction analysis.

    All query and transformation methods follow immutable patterns — the
    internally stored DataFrame is never modified after loading.

    Attributes:
        config: Optional configuration dictionary for runtime settings.
        _df: Internal DataFrame holding the loaded interaction data,
             or None if no data has been loaded yet.

    Example:
        >>> db = DrugInteractionDB()
        >>> db.load_data("demo/sample_data.csv")
        >>> results = db.lookup("warfarin", "aspirin")
        >>> print(results[["drug_a", "drug_b", "severity", "clinical_effect"]])
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None) -> None:
        """
        Initialise the DrugInteractionDB instance.

        Args:
            config: Optional dictionary with runtime configuration options.
                    Supported keys:
                      - ``strict_validation`` (bool): Raise on unknown columns.
        """
        self.config: Dict[str, Any] = config if config is not None else {}
        self._df: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------

    def load_data(self, filepath: str) -> pd.DataFrame:
        """
        Load interaction data from a CSV or Excel file.

        The loaded DataFrame is also stored internally so subsequent query
        methods can be called without passing it explicitly.

        Args:
            filepath: Absolute or relative path to a ``.csv``, ``.xlsx``,
                      or ``.xls`` file.

        Returns:
            Preprocessed DataFrame with standardised column names.

        Raises:
            FileNotFoundError: If *filepath* does not exist.
            ValueError: If the file is empty or lacks required columns.
            TypeError: If *filepath* is not a string.
        """
        if not isinstance(filepath, str):
            raise TypeError(
                f"'filepath' must be a string, got {type(filepath).__name__!r}."
            )

        p = Path(filepath)
        if not p.exists():
            raise FileNotFoundError(f"Data file not found: {filepath}")

        if p.suffix in (".xlsx", ".xls"):
            raw: pd.DataFrame = pd.read_excel(filepath)
        else:
            raw = pd.read_csv(filepath)

        processed = self.preprocess(raw)
        self.validate(processed)
        self._df = processed
        return processed

    # ------------------------------------------------------------------
    # Validation & preprocessing
    # ------------------------------------------------------------------

    def validate(self, df: pd.DataFrame) -> bool:
        """
        Validate that *df* satisfies minimum structural requirements.

        Checks that the DataFrame is non-empty and contains all columns
        listed in :data:`REQUIRED_COLUMNS`.

        Args:
            df: DataFrame to validate.

        Returns:
            ``True`` when validation passes.

        Raises:
            TypeError: If *df* is not a pandas DataFrame.
            ValueError: If *df* is empty or missing required columns.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f"Expected a pandas DataFrame, got {type(df).__name__!r}."
            )
        if df.empty:
            raise ValueError("Input DataFrame is empty.")

        missing = REQUIRED_COLUMNS - set(df.columns)
        if missing:
            raise ValueError(
                f"Missing required columns: {sorted(missing)}. "
                f"Present columns: {sorted(df.columns)}"
            )
        return True

    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean and normalise raw interaction data.

        All operations are immutable — the original DataFrame is never mutated.
        Steps performed:

        1. Drop fully-empty rows.
        2. Standardise column names to lowercase with underscores.
        3. Strip leading/trailing whitespace from string columns.
        4. Normalise ``drug_a``, ``drug_b``, ``severity``, and
           ``interaction_type`` to lowercase for case-insensitive matching.
        5. Drop exact duplicate rows (bidirectional duplicates where
           ``(A, B)`` and ``(B, A)`` refer to the same pair are retained as
           distinct records since directionality may carry meaning in source
           data).

        Args:
            df: Raw DataFrame as loaded from the source file.

        Returns:
            New cleaned DataFrame; the input is never modified.

        Raises:
            TypeError: If *df* is not a pandas DataFrame.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f"Expected a pandas DataFrame, got {type(df).__name__!r}."
            )

        result = df.copy()
        result = result.dropna(how="all")
        result = result.rename(
            columns={c: c.lower().strip().replace(" ", "_") for c in result.columns}
        )

        # Strip whitespace from all string columns
        str_cols: List[str] = list(result.select_dtypes(include="object").columns)
        result = result.assign(
            **{col: result[col].str.strip() for col in str_cols}
        )

        # Normalise key identifier columns to lowercase for case-insensitive matching
        key_cols: tuple = ("drug_a", "drug_b", "severity", "interaction_type")
        for col in key_cols:
            if col in result.columns:
                result = result.assign(**{col: result[col].str.lower()})

        # Remove exact duplicates
        result = result.drop_duplicates()

        return result

    # ------------------------------------------------------------------
    # Query methods
    # ------------------------------------------------------------------

    def _require_loaded(self) -> pd.DataFrame:
        """
        Return the internally stored DataFrame, raising if not yet loaded.

        Returns:
            The preprocessed DataFrame stored by :meth:`load_data`.

        Raises:
            RuntimeError: If :meth:`load_data` has not been called yet.
        """
        if self._df is None:
            raise RuntimeError(
                "No data loaded. Call load_data(filepath) first."
            )
        return self._df

    def lookup(
        self,
        drug_a: str,
        drug_b: str,
        df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Look up all interactions between *drug_a* and *drug_b*.

        The search is bidirectional and case-insensitive: an interaction
        recorded as (A, B) will also be returned when querying (B, A).
        Leading/trailing whitespace in drug names is ignored.

        Args:
            drug_a: Name of the first drug. Must be a non-empty string.
            drug_b: Name of the second drug. Must be a non-empty string.
            df: Optional DataFrame to query. When omitted the internally
                loaded DataFrame is used.

        Returns:
            DataFrame of matching interaction rows (may be empty if no
            interaction is recorded for this pair).

        Raises:
            TypeError: If either drug name is not a string.
            ValueError: If either drug name is empty or whitespace-only.
            RuntimeError: If no DataFrame is available (neither *df* nor a
                          previously loaded one).
        """
        a = _validate_drug_name(drug_a, "drug_a")
        b = _validate_drug_name(drug_b, "drug_b")
        source = df if df is not None else self._require_loaded()

        mask = (
            ((source["drug_a"] == a) & (source["drug_b"] == b))
            | ((source["drug_a"] == b) & (source["drug_b"] == a))
        )
        return source.loc[mask].copy()

    def filter_by_severity(
        self,
        severity: str,
        df: Optional[pd.DataFrame] = None,
        minimum: bool = False,
    ) -> pd.DataFrame:
        """
        Filter interactions by severity level.

        Args:
            severity: Target severity — one of ``minor``, ``moderate``,
                      ``major``, or ``contraindicated`` (case-insensitive).
            df: Optional DataFrame to filter. Falls back to the loaded data.
            minimum: When ``True``, return all interactions at *severity* or
                     above in the clinical risk scale (e.g. passing
                     ``"major"`` returns both ``major`` and
                     ``contraindicated`` rows).

        Returns:
            Filtered DataFrame (a copy; the source is not mutated).

        Raises:
            TypeError: If *severity* is not a string.
            ValueError: If *severity* is not a recognised level.
            RuntimeError: If no data is available.
        """
        if not isinstance(severity, str):
            raise TypeError(
                f"'severity' must be a string, got {type(severity).__name__!r}."
            )
        source = df if df is not None else self._require_loaded()
        level = severity.lower().strip()

        if level not in SEVERITY_ORDER:
            raise ValueError(
                f"Unknown severity '{severity}'. "
                f"Choose from: {SEVERITY_ORDER}"
            )

        if minimum:
            threshold = SEVERITY_ORDER.index(level)
            allowed = frozenset(SEVERITY_ORDER[threshold:])
            return source.loc[source["severity"].isin(allowed)].copy()

        return source.loc[source["severity"] == level].copy()

    def search_drug(
        self,
        drug: str,
        df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Return all interactions that involve *drug* in either position.

        The search is case-insensitive and strips surrounding whitespace.

        Args:
            drug: Drug name to search for. Must be a non-empty string.
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            DataFrame containing all rows where *drug* appears as
            ``drug_a`` or ``drug_b``.

        Raises:
            TypeError: If *drug* is not a string.
            ValueError: If *drug* is empty or whitespace-only.
            RuntimeError: If no data is available.
        """
        name = _validate_drug_name(drug, "drug")
        source = df if df is not None else self._require_loaded()
        mask = (source["drug_a"] == name) | (source["drug_b"] == name)
        return source.loc[mask].copy()

    def check_multi_drug(
        self,
        drugs: List[str],
        df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Check all pairwise interactions within a list of drugs.

        For *n* drugs this produces up to *n(n-1)/2* unique pairs, each
        resolved via :meth:`lookup`. Duplicate result rows (same interaction
        triggered by different queries) are deduplicated in the output.

        Args:
            drugs: List of drug names to check (case-insensitive). Must
                   contain at least two entries. Each entry must be a
                   non-empty string.
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            DataFrame of all found interactions with an extra
            ``queried_pair`` column showing which pair triggered each row.
            Returns an empty DataFrame if no interactions are found.

        Raises:
            TypeError: If *drugs* is not a list or contains non-string items.
            ValueError: If fewer than two drugs are provided or any drug
                        name is empty.
            RuntimeError: If no data is available.
        """
        if not isinstance(drugs, list):
            raise TypeError(
                f"'drugs' must be a list, got {type(drugs).__name__!r}."
            )
        if len(drugs) < 2:
            raise ValueError("Provide at least two drugs for multi-drug check.")

        source = df if df is not None else self._require_loaded()
        frames: List[pd.DataFrame] = []

        for i, a in enumerate(drugs):
            a_validated = _validate_drug_name(a, f"drugs[{i}]")
            for j, b in enumerate(drugs[i + 1:], start=i + 1):
                b_validated = _validate_drug_name(b, f"drugs[{j}]")
                hits = self.lookup(a_validated, b_validated, df=source)
                if not hits.empty:
                    labeled = hits.assign(
                        queried_pair=f"{a_validated} / {b_validated}"
                    )
                    frames.append(labeled)

        if not frames:
            return pd.DataFrame()

        combined = pd.concat(frames, ignore_index=True)
        dedup_cols = [c for c in combined.columns if c != "queried_pair"]
        return combined.drop_duplicates(subset=dedup_cols)

    # ------------------------------------------------------------------
    # Analysis helpers
    # ------------------------------------------------------------------

    def analyze(self, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Run summary analysis on an interaction DataFrame.

        Computes record counts, missing-value percentages, severity
        distribution, and interaction-type distribution. The input
        DataFrame is preprocessed before analysis but never mutated.

        Args:
            df: Preprocessed or raw interaction DataFrame.

        Returns:
            Dictionary with the following keys:

            - ``total_records`` (int): Number of rows after preprocessing.
            - ``columns`` (list of str): Column names present.
            - ``missing_pct`` (dict str → float): Percentage of missing
              values per column, rounded to one decimal place.
            - ``severity_counts`` (dict str → int): Row count per severity
              level (only present when the ``severity`` column exists).
            - ``interaction_type_counts`` (dict str → int): Row count per
              interaction type (only present when ``interaction_type``
              column exists).
            - ``summary_stats`` (dict): Descriptive statistics for numeric
              columns (only present when numeric columns exist).

        Raises:
            TypeError: If *df* is not a pandas DataFrame.
        """
        cleaned = self.preprocess(df)
        result: Dict[str, Any] = {
            "total_records": len(cleaned),
            "columns": list(cleaned.columns),
            "missing_pct": (
                cleaned.isnull().sum() / max(len(cleaned), 1) * 100
            ).round(1).to_dict(),
        }

        if "severity" in cleaned.columns:
            result["severity_counts"] = (
                cleaned["severity"].value_counts().to_dict()
            )

        if "interaction_type" in cleaned.columns:
            result["interaction_type_counts"] = (
                cleaned["interaction_type"].value_counts().to_dict()
            )

        numeric_df = cleaned.select_dtypes(include="number")
        if not numeric_df.empty:
            result["summary_stats"] = numeric_df.describe().round(3).to_dict()

        return result

    def run(self, filepath: str) -> Dict[str, Any]:
        """
        Convenience pipeline: load a file, validate it, and return analysis.

        Equivalent to calling :meth:`load_data` followed by :meth:`analyze`.

        Args:
            filepath: Path to the CSV or Excel interaction data file.

        Returns:
            Analysis result dictionary from :meth:`analyze`.

        Raises:
            FileNotFoundError: If *filepath* does not exist.
            TypeError: If *filepath* is not a string.
            ValueError: If the data is empty or missing required columns.
        """
        df = self.load_data(filepath)
        return self.analyze(df)

    def to_dataframe(self, result: Dict[str, Any]) -> pd.DataFrame:
        """
        Flatten an analysis result dictionary into a two-column DataFrame.

        Nested dictionaries are expanded with dotted keys, e.g.
        ``severity_counts.major``. Non-dict scalar values are kept as-is.

        Args:
            result: Dictionary as returned by :meth:`analyze`.

        Returns:
            DataFrame with columns ``metric`` and ``value``.

        Raises:
            TypeError: If *result* is not a dictionary.
        """
        if not isinstance(result, dict):
            raise TypeError(
                f"'result' must be a dict, got {type(result).__name__!r}."
            )
        rows: List[Dict[str, Any]] = []
        for k, v in result.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    rows.append({"metric": f"{k}.{kk}", "value": vv})
            else:
                rows.append({"metric": k, "value": v})
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    def list_drugs(self, df: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Return a sorted list of all unique drug names in the dataset.

        Drug names are collected from both the ``drug_a`` and ``drug_b``
        columns, deduplicated, and returned in ascending alphabetical order.

        Args:
            df: Optional DataFrame to inspect. Falls back to the loaded data.

        Returns:
            Sorted list of lowercase drug names appearing in either column.

        Raises:
            RuntimeError: If no data is available.
        """
        source = df if df is not None else self._require_loaded()
        drugs = pd.concat(
            [source["drug_a"], source["drug_b"]], ignore_index=True
        ).dropna().unique()
        return sorted(str(d) for d in drugs)

    def get_contraindicated(
        self, df: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Return all interactions flagged as contraindicated.

        Convenience wrapper around :meth:`filter_by_severity` for the
        highest-risk category.

        Args:
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            DataFrame of contraindicated interaction rows.

        Raises:
            RuntimeError: If no data is available.
        """
        return self.filter_by_severity("contraindicated", df=df)

    def has_interaction(
        self,
        drug_a: str,
        drug_b: str,
        df: Optional[pd.DataFrame] = None,
    ) -> bool:
        """
        Return True if any interaction exists between *drug_a* and *drug_b*.

        Convenience wrapper around :meth:`lookup` for boolean checks.

        Args:
            drug_a: Name of the first drug. Must be a non-empty string.
            drug_b: Name of the second drug. Must be a non-empty string.
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            ``True`` if at least one interaction row is found, else ``False``.

        Raises:
            TypeError: If either drug name is not a string.
            ValueError: If either drug name is empty or whitespace-only.
            RuntimeError: If no data is available.
        """
        return not self.lookup(drug_a, drug_b, df=df).empty

    def get_high_risk_drugs(
        self, df: Optional[pd.DataFrame] = None
    ) -> List[str]:
        """
        Return a sorted list of drugs that appear in major or contraindicated interactions.

        Useful for quickly identifying drugs that require heightened clinical
        vigilance in prescribing.

        Args:
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            Sorted list of lowercase drug names involved in interactions
            with severity ``major`` or ``contraindicated``.

        Raises:
            RuntimeError: If no data is available.
        """
        high_risk = self.filter_by_severity("major", df=df, minimum=True)
        if high_risk.empty:
            return []
        drugs = pd.concat(
            [high_risk["drug_a"], high_risk["drug_b"]], ignore_index=True
        ).dropna().unique()
        return sorted(str(d) for d in drugs)
