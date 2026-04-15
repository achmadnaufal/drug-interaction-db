# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
