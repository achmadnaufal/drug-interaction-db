# Drug Interaction DB

Drug-drug interaction database query and risk assessment tools.

## Features

- Load interaction data from CSV or Excel files
- Case-insensitive, bidirectional drug lookup
- Severity filtering with threshold support (`minor` → `contraindicated`)
- Multi-drug pairwise interaction check
- Summary analysis with severity and type distributions
- Immutable data transformations — source data is never mutated
- Comprehensive input validation with clear error messages

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

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
```

## Sample Data

The file `demo/sample_data.csv` contains 20 clinically representative
drug-drug interactions. Columns:

| Column | Description |
|---|---|
| `drug_a` | First drug name |
| `drug_b` | Second drug name |
| `interaction_type` | `pharmacokinetic` or `pharmacodynamic` |
| `severity` | `minor`, `moderate`, `major`, or `contraindicated` |
| `mechanism` | Pharmacological mechanism |
| `clinical_effect` | Expected patient outcome |
| `recommendation` | Clinical management advice |
| `evidence_level` | `A` (strong) or `B` (moderate) |
| `source` | Literature citation |

Example rows:

```
drug_a,drug_b,interaction_type,severity,...
warfarin,aspirin,pharmacodynamic,major,...
ssris,maois,pharmacodynamic,contraindicated,...
sildenafil,nitrates,pharmacodynamic,contraindicated,...
```

## Example Queries

```python
# Find all interactions involving amiodarone
amio = db.search_drug("amiodarone")

# List all drugs in the database
all_drugs = db.list_drugs()

# Get all contraindicated pairs
danger = db.get_contraindicated()

# Full analysis summary
summary = db.analyze(db._df)
print(summary["severity_counts"])

# Run entire pipeline from file path
result = db.run("demo/sample_data.csv")
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
tests/test_interactions.py::TestLookup::test_lookup_known_pair_returns_rows PASSED
tests/test_interactions.py::TestLookup::test_lookup_known_pair_correct_severity PASSED
...
PASSED (30+ tests)
```

## Data Format

Expected CSV columns (minimum required): `drug_a`, `drug_b`, `interaction_type`, `severity`

All columns are normalised to lowercase on load. Drug names are matched
case-insensitively and leading/trailing whitespace is stripped automatically.

## Project Structure

```
drug-interaction-db/
├── src/
│   ├── __init__.py
│   ├── main.py           # Core DrugInteractionDB class
│   └── data_generator.py # Synthetic data generator
├── demo/
│   └── sample_data.csv   # 20-row realistic interaction dataset
├── tests/
│   ├── __init__.py
│   └── test_interactions.py  # pytest unit tests
├── examples/
│   └── basic_usage.py    # End-to-end usage example
├── data/                 # Drop real datasets here (gitignored)
├── CHANGELOG.md
├── requirements.txt
└── README.md
```

## License

MIT License — free to use, modify, and distribute.
