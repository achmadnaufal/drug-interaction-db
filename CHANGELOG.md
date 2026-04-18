# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased] - 2026-04-19

### Added
- `src/alternative_suggester.py`: new module that proposes safer replacement
  drugs from the same therapeutic class when a patient's regimen contains a
  drug causing a problematic interaction.  Ships a curated
  `THERAPEUTIC_CLASSES` registry (PPIs, SSRIs, statins, macrolides, azole
  antifungals, NSAIDs, ACE inhibitors, anticoagulants) and scores each
  candidate against the rest of the regimen using the severity weights from
  the polypharmacy scorer.
- `suggest_alternatives(target, regimen, db, ...)`: returns an immutable
  tuple of `AlternativeSuggestion` objects sorted by ascending DDI burden.
  Handles unknown drugs (returns empty tuple), empty regimens (raises
  `ValueError`), duplicate entries (deduplicated internally), and
  case-insensitive matching.
- `suggestions_to_dataframe(suggestions)`: flattens a suggestions tuple into
  a pandas DataFrame for reporting.
- `list_classes()`: exposes the registered therapeutic class names.
- `AlternativeSuggestion` frozen dataclass with `candidate`,
  `therapeutic_class`, `total_score`, `highest_severity`,
  `interacting_drugs`, and `is_safer` fields.
- `tests/test_alternative_suggester.py`: 27 pytest tests covering happy
  path, case-insensitivity, unknown drugs, empty / duplicate regimens,
  max-results capping, `is_safer` flag semantics, full input validation,
  DataFrame round-tripping, and class-registry invariants.
- `src/__init__.py`: re-exports the new public API.
- README top-level "Disclaimer — Not Medical Advice" block, a
  "Severity Classification" table, and a "New: Alternative-Drug Suggester"
  section with step-by-step usage, example output, and an API reference.

## [Unreleased] - 2026-04-18

### Added
- `src/cyp450_tagger.py`: new module that scans the `mechanism` text of every
  DDI record and tags each row with the CYP450 isoenzymes (and related
  transporters such as P-glycoprotein) implicated in the interaction.
- `extract_enzymes(mechanism)`: pure function returning a tuple of immutable
  `EnzymeTag` objects, with case-insensitive, word-boundary-aware matching.
- `tag_interactions(df)`: appends `cyp_enzymes` and `cyp_role` columns to a
  DDI DataFrame without mutating the input.
- `summarise_by_enzyme(df)`: aggregates inhibitor / inducer / unspecified
  counts per enzyme, sorted by total descending.
- `filter_by_enzyme(db, enzyme, role=None)`: returns interaction rows that
  mention a given CYP isoform, optionally restricted by metabolic role.
- `EnzymeTag` frozen dataclass and `KNOWN_ENZYMES` constant tuple covering
  common CYP isoforms (CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, ...) plus
  P-glycoprotein, UGT1A1, and OATP1B1.
- `tests/test_cyp450_tagger.py`: 36 pytest tests covering happy path,
  empty/whitespace/None input, NaN handling, missing column, unknown
  enzyme/role, case sensitivity, word-boundary protection, role precedence
  (inhibition over induction), and immutability of input DataFrames.
- `src/__init__.py`: re-exports the new public API alongside the existing
  `DrugInteractionDB` and polypharmacy scorer entry points.
- README "New: CYP450 Enzyme-Mediated Interaction Tagger" section with
  step-by-step usage, sample output, and API reference table.

## [Unreleased] - 2026-04-17

### Added
- `src/polypharmacy_risk_scorer.py`: new module that scores a patient medication
  list by total DDI burden — enumerates all pairwise interactions, applies
  severity weights (`minor=1`, `moderate=3`, `major=7`, `contraindicated=15`),
  sums into a `total_score`, and classifies regimens as `"low"`, `"high"`, or
  `"critical"`.  Flags polypharmacy when 5 or more drugs are present.
- `RiskReport` and `PairScore` frozen dataclasses for immutable, type-safe
  result objects.
- `score_regimen()` function with full input validation, deduplication, and
  configurable thresholds and weights.
- `report_to_dataframe()` helper to flatten a `RiskReport` into a sorted
  pair-score DataFrame.
- `tests/test_polypharmacy_risk_scorer.py`: 29 pytest tests covering happy path,
  single-drug, unknown drugs, duplicate deduplication, contraindicated flagging,
  polypharmacy threshold, input validation, determinism, and custom weights.
- README "New: Polypharmacy Risk Scorer" section with step-by-step usage,
  sample output, custom threshold examples, and API reference table.

## [0.2.0] - 2026-04-16

### Added
- `has_interaction(drug_a, drug_b)` convenience method returning a boolean
- `get_high_risk_drugs()` method listing drugs involved in major or contraindicated interactions
- Full type hints on every function and method signature
- Comprehensive docstrings covering Args, Returns, and Raises for every public method
- `_validate_drug_name()` module-level helper for consistent input validation
- Input validation raises `TypeError` for non-string drug names and `ValueError` for empty or whitespace-only names
- `TypeError` guards on `validate()`, `preprocess()`, `filter_by_severity()`, `to_dataframe()`, `check_multi_drug()`, and `load_data()`
- `REQUIRED_COLUMNS` changed from `set` to `frozenset` to enforce immutability
- `frozenset` used for severity threshold set in `filter_by_severity()` (immutable pattern)
- Expanded test suite: 60+ pytest unit tests across 10 test classes covering all public methods, edge cases, and type error paths
- `TestConstants`, `TestHasInteraction`, `TestGetHighRiskDrugs`, `TestLoadData`, and `TestInit` test classes added
- `data_generator.py` rewritten with curated pharmacologically accurate interactions, type hints, and full docstrings
- README updated with badges, Sample Output section, API Reference table, and improved Project Structure

### Changed
- `data_generator.py` now produces realistic drug-drug interaction records instead of random numeric data
- `preprocess()` uses `rename(columns=...)` instead of direct column assignment for stricter immutability
- `analyze()` guards against division-by-zero on empty DataFrames via `max(len(cleaned), 1)`
- README restructured with Quick Start, Sample Output, and full API Reference sections

### Fixed
- `check_multi_drug()` now validates each individual drug name in the list, not just the list length
- `preprocess()` raises `TypeError` clearly when given a non-DataFrame argument

## [0.1.0] - 2024-01-01

### Added
- Initial project scaffold
- `DrugInteractionDB` core class with `load_data`, `validate`, `preprocess`, and `analyze` methods
- `lookup`, `filter_by_severity`, `search_drug`, `check_multi_drug` query methods
- `list_drugs` and `get_contraindicated` utility helpers
- Sample data generator script
- Basic usage example
