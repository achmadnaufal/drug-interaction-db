# Drug Interaction DB

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Tests](https://img.shields.io/badge/tests-pytest-orange)
![Coverage](https://img.shields.io/badge/coverage-80%25%2B-brightgreen)
![Code Style](https://img.shields.io/badge/code%20style-immutable-lightgrey)

Drug-drug interaction database query and risk assessment tools.

## Features

- Load interaction data from CSV or Excel files
- Case-insensitive, bidirectional drug lookup
- Severity filtering with threshold support (`minor` to `contraindicated`)
- Multi-drug pairwise interaction check across a full medication list
- Summary analysis with severity and type distributions
- Immutable data transformations — source data is never mutated
- Comprehensive input validation with clear error messages
- Type-safe API with full type hints throughout

## Quick Start

```bash
pip install -r requirements.txt
```

```python
from src.main import DrugInteractionDB

db = DrugInteractionDB()
db.load_data("demo/sample_data.csv")

# Look up a specific drug pair (bidirectional, case-insensitive)
hits = db.lookup("Warfarin", "Aspirin")
print(hits[["drug_a", "drug_b", "severity", "clinical_effect"]])

# Filter all major interactions
major = db.filter_by_severity("major")
print(f"Major interactions: {len(major)}")

# All interactions at 'major' severity or above
high_risk = db.filter_by_severity("major", minimum=True)

# Check a patient's full medication list for pairwise interactions
result = db.check_multi_drug(["warfarin", "aspirin", "omeprazole"])
print(result[["queried_pair", "severity", "recommendation"]])

# Boolean check
if db.has_interaction("sildenafil", "nitrates"):
    print("Contraindicated combination detected!")

# Drugs requiring heightened clinical vigilance
risky = db.get_high_risk_drugs()
print(f"High-risk drugs in dataset: {risky}")
```

## Sample Output

```
   drug_a  drug_b severity                             clinical_effect
  warfarin aspirin    major  Increased risk of bleeding including ...

Major interactions: 13

          queried_pair severity                          recommendation
  warfarin / aspirin    major  Avoid combination; if unavoidable ...

Contraindicated combination detected!
High-risk drugs in dataset: ['ace_inhibitors', 'allopurinol', 'amiodarone', ...]
```

## Sample Data

The file `demo/sample_data.csv` contains 20 clinically representative
drug-drug interactions drawn from published pharmacology literature.

| Column | Description |
|---|---|
| `drug_a` | First drug name |
| `drug_b` | Second drug name |
| `interaction_type` | `pharmacokinetic` or `pharmacodynamic` |
| `severity` | `minor`, `moderate`, `major`, or `contraindicated` |
| `mechanism` | Pharmacological mechanism of interaction |
| `clinical_effect` | Expected patient outcome |
| `recommendation` | Clinical management advice |
| `evidence_level` | `A` (strong) or `B` (moderate) |
| `source` | Literature citation |

Selected interactions included:

| Drug A | Drug B | Severity | Mechanism |
|---|---|---|---|
| Warfarin | Aspirin | major | Additive anticoagulant + antiplatelet |
| SSRIs | MAOIs | contraindicated | Serotonin syndrome |
| Sildenafil | Nitrates | contraindicated | Additive vasodilation |
| Simvastatin | Amiodarone | major | CYP3A4 inhibition |
| Clopidogrel | Omeprazole | moderate | CYP2C19 inhibition |
| Allopurinol | Azathioprine | major | Xanthine oxidase inhibition |
| Fluoxetine | Tamoxifen | major | CYP2D6 inhibition |

## Example Queries

```python
# Find all interactions involving amiodarone
amio = db.search_drug("amiodarone")

# List all drugs in the database
all_drugs = db.list_drugs()

# Get all contraindicated pairs
danger = db.get_contraindicated()
print(f"Contraindicated pairs: {len(danger)}")

# Full analysis summary
summary = db.analyze(db._df)
print(summary["severity_counts"])
# {'major': 13, 'moderate': 4, 'contraindicated': 2, 'minor': 1}

# Flatten analysis to a tabular DataFrame
flat = db.to_dataframe(summary)

# Run entire pipeline from file path
result = db.run("demo/sample_data.csv")
```

## Data Format

Minimum required CSV columns: `drug_a`, `drug_b`, `interaction_type`, `severity`

All columns are normalised to lowercase on load. Drug names are matched
case-insensitively and leading/trailing whitespace is stripped automatically.
Exact duplicate rows are removed during preprocessing.

## Project Structure

```
drug-interaction-db/
├── src/
│   ├── __init__.py
│   ├── main.py              # Core DrugInteractionDB class
│   └── data_generator.py    # Curated sample data generator
├── demo/
│   └── sample_data.csv      # 20-row realistic interaction dataset
├── tests/
│   ├── __init__.py
│   └── test_interactions.py # pytest unit tests (80%+ coverage)
├── examples/
│   └── basic_usage.py       # End-to-end usage example
├── data/                    # Drop real datasets here (gitignored)
├── CHANGELOG.md
├── requirements.txt
└── README.md
```

## Running Tests

```bash
# Install dependencies
pip install -r requirements.txt

# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ -v --cov=src --cov-report=term-missing
```

Expected output:

```
tests/test_interactions.py::TestConstants::test_severity_order_has_four_levels PASSED
tests/test_interactions.py::TestLookup::test_lookup_known_pair_returns_rows PASSED
tests/test_interactions.py::TestLookup::test_lookup_known_pair_correct_severity PASSED
tests/test_interactions.py::TestSeverityFilter::test_filter_major_returns_only_major PASSED
...
========= 60+ passed in 0.5s =========
```

## API Reference

| Method | Description |
|---|---|
| `load_data(filepath)` | Load CSV or Excel file and store internally |
| `lookup(drug_a, drug_b)` | Bidirectional search for a specific pair |
| `has_interaction(drug_a, drug_b)` | Boolean check for any interaction |
| `filter_by_severity(severity, minimum=False)` | Filter by severity level |
| `search_drug(drug)` | All interactions involving a single drug |
| `check_multi_drug(drugs)` | Pairwise check across a medication list |
| `get_contraindicated()` | All contraindicated pairs |
| `get_high_risk_drugs()` | Drugs in major or contraindicated pairs |
| `list_drugs()` | Sorted list of all unique drug names |
| `analyze(df)` | Summary statistics dictionary |
| `to_dataframe(result)` | Flatten analysis dict to two-column DataFrame |
| `run(filepath)` | Load + analyze pipeline in one call |

## License

MIT License — free to use, modify, and distribute.
