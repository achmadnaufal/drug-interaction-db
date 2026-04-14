"""
Drug-drug interaction database query and risk assessment tools.

This module provides the DrugInteractionDB class for loading, validating,
querying, and analyzing drug-drug interaction data from CSV or Excel files.

Author: github.com/achmadnaufal
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple


SEVERITY_ORDER = ["minor", "moderate", "major", "contraindicated"]

REQUIRED_COLUMNS = {"drug_a", "drug_b", "interaction_type", "severity"}


class DrugInteractionDB:
    """
    Drug interaction database query and risk assessment tool.

    Loads drug-drug interaction data from CSV or Excel files and provides
    methods for lookup, severity filtering, bidirectional search, and
    multi-drug interaction analysis.

    Attributes:
        config: Optional configuration dictionary for runtime settings.
        _df: Internal DataFrame holding the loaded interaction data,
             or None if no data has been loaded yet.

    Example:
        >>> db = DrugInteractionDB()
        >>> db.load_data("demo/sample_data.csv")
        >>> results = db.lookup("warfarin", "aspirin")
    """

    def __init__(self, config: Optional[Dict] = None) -> None:
        """
        Initialise the DrugInteractionDB instance.

        Args:
            config: Optional dictionary with runtime configuration options.
                    Supported keys:
                      - ``strict_validation`` (bool): Raise on unknown columns.
        """
        self.config: Dict = config or {}
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
        """
        p = Path(filepath)
        if not p.exists():
            raise FileNotFoundError(f"Data file not found: {filepath}")

        if p.suffix in (".xlsx", ".xls"):
            raw = pd.read_excel(filepath)
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

        Args:
            df: DataFrame to validate.

        Returns:
            ``True`` when validation passes.

        Raises:
            ValueError: If *df* is empty or missing required columns.
        """
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

        Operations performed (all immutable — original DataFrame is not mutated):
          - Drop fully-empty rows.
          - Standardise column names to lowercase with underscores.
          - Strip leading/trailing whitespace from string columns.
          - Normalise ``drug_a``, ``drug_b``, and ``severity`` to lowercase.
          - Drop exact duplicate rows.

        Args:
            df: Raw DataFrame as loaded from the source file.

        Returns:
            New cleaned DataFrame; the input is never modified.
        """
        result = df.copy()
        result = result.dropna(how="all")
        result.columns = [c.lower().strip().replace(" ", "_") for c in result.columns]

        # Strip whitespace from all string columns
        str_cols = result.select_dtypes(include="object").columns
        result = result.assign(
            **{col: result[col].str.strip() for col in str_cols}
        )

        # Normalise key identifier columns to lowercase for case-insensitive matching
        for col in ("drug_a", "drug_b", "severity", "interaction_type"):
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
            The preprocessed DataFrame stored by ``load_data``.

        Raises:
            RuntimeError: If ``load_data`` has not been called yet.
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

        Args:
            drug_a: Name of the first drug.
            drug_b: Name of the second drug.
            df: Optional DataFrame to query. When omitted the internally
                loaded DataFrame is used.

        Returns:
            DataFrame of matching interaction rows (may be empty).

        Raises:
            RuntimeError: If no DataFrame is available (neither *df* nor a
                          previously loaded one).
        """
        source = df if df is not None else self._require_loaded()
        a = drug_a.lower().strip()
        b = drug_b.lower().strip()

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
                     above in the clinical risk scale.

        Returns:
            Filtered DataFrame.

        Raises:
            ValueError: If *severity* is not a recognised level.
            RuntimeError: If no data is available.
        """
        source = df if df is not None else self._require_loaded()
        level = severity.lower().strip()

        if level not in SEVERITY_ORDER:
            raise ValueError(
                f"Unknown severity '{severity}'. "
                f"Choose from: {SEVERITY_ORDER}"
            )

        if minimum:
            threshold = SEVERITY_ORDER.index(level)
            allowed = set(SEVERITY_ORDER[threshold:])
            return source.loc[source["severity"].isin(allowed)].copy()

        return source.loc[source["severity"] == level].copy()

    def search_drug(
        self,
        drug: str,
        df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Return all interactions that involve *drug* in either position.

        Args:
            drug: Drug name to search for (case-insensitive).
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            DataFrame containing all rows where *drug* appears as
            ``drug_a`` or ``drug_b``.

        Raises:
            RuntimeError: If no data is available.
        """
        source = df if df is not None else self._require_loaded()
        name = drug.lower().strip()
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
        resolved via :meth:`lookup`.

        Args:
            drugs: List of drug names to check (case-insensitive).
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            DataFrame of all found interactions, deduplicated, with an extra
            ``queried_pair`` column showing which pair triggered each row.

        Raises:
            ValueError: If fewer than two drugs are provided.
            RuntimeError: If no data is available.
        """
        if len(drugs) < 2:
            raise ValueError("Provide at least two drugs for multi-drug check.")

        source = df if df is not None else self._require_loaded()
        frames: List[pd.DataFrame] = []

        for i, a in enumerate(drugs):
            for b in drugs[i + 1 :]:
                hits = self.lookup(a, b, df=source)
                if not hits.empty:
                    labeled = hits.assign(queried_pair=f"{a.lower()} / {b.lower()}")
                    frames.append(labeled)

        if not frames:
            return pd.DataFrame()

        combined = pd.concat(frames, ignore_index=True)
        return combined.drop_duplicates(subset=[c for c in combined.columns if c != "queried_pair"])

    # ------------------------------------------------------------------
    # Analysis helpers
    # ------------------------------------------------------------------

    def analyze(self, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Run summary analysis on an interaction DataFrame.

        Computes record counts, missing-value percentages, severity
        distribution, and interaction-type distribution.

        Args:
            df: Preprocessed interaction DataFrame.

        Returns:
            Dictionary with keys:
              - ``total_records`` (int)
              - ``columns`` (list of str)
              - ``missing_pct`` (dict col → float)
              - ``severity_counts`` (dict severity → int, if column present)
              - ``interaction_type_counts`` (dict type → int, if column present)
              - ``summary_stats`` (dict, numeric columns only)
        """
        cleaned = self.preprocess(df)
        result: Dict[str, Any] = {
            "total_records": len(cleaned),
            "columns": list(cleaned.columns),
            "missing_pct": (
                cleaned.isnull().sum() / len(cleaned) * 100
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

        Args:
            filepath: Path to the CSV or Excel interaction data file.

        Returns:
            Analysis result dictionary from :meth:`analyze`.

        Raises:
            FileNotFoundError: If *filepath* does not exist.
            ValueError: If the data is empty or missing required columns.
        """
        df = self.load_data(filepath)
        return self.analyze(df)

    def to_dataframe(self, result: Dict) -> pd.DataFrame:
        """
        Flatten an analysis result dictionary into a two-column DataFrame.

        Nested dictionaries are expanded with dotted keys, e.g.
        ``severity_counts.major``.

        Args:
            result: Dictionary as returned by :meth:`analyze`.

        Returns:
            DataFrame with columns ``metric`` and ``value``.
        """
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

        Args:
            df: Optional DataFrame to inspect. Falls back to the loaded data.

        Returns:
            Sorted list of lowercase drug names appearing in either
            ``drug_a`` or ``drug_b`` columns.

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

        Args:
            df: Optional DataFrame to query. Falls back to the loaded data.

        Returns:
            DataFrame of contraindicated interaction rows.

        Raises:
            RuntimeError: If no data is available.
        """
        return self.filter_by_severity("contraindicated", df=df)
