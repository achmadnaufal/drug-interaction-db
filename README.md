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

## New: Polypharmacy Risk Scorer

The `src/polypharmacy_risk_scorer.py` module quantifies the total DDI burden for
a patient's full medication list by scoring every pairwise combination and
aggregating into a single risk report.

### How It Works

1. Every unique drug pair is looked up in the loaded interaction database.
2. The highest severity level for each pair is mapped to a numeric weight
   (`minor=1`, `moderate=3`, `major=7`, `contraindicated=15`).
3. Pair scores are summed into a `total_score`.
4. The regimen is classified as `"low"`, `"high"` (score ≥ 10), or
   `"critical"` (score ≥ 25 or any contraindicated pair).
5. A `polypharmacy_flagged` flag is set when 5 or more drugs are present.

### Step-by-Step Usage

```python
from src.main import DrugInteractionDB
from src.polypharmacy_risk_scorer import score_regimen, report_to_dataframe

# 1. Load the interaction database
db = DrugInteractionDB()
db.load_data("demo/sample_data.csv")

# 2. Define the patient's medication list
patient_meds = ["warfarin", "aspirin", "fluconazole", "simvastatin", "amiodarone"]

# 3. Score the regimen
report = score_regimen(patient_meds, db)

# 4. Inspect the summary
print(f"Total DDI burden score : {report.total_score}")
print(f"Risk level             : {report.risk_level}")
print(f"Interacting pairs      : {report.interacting_pairs} / {report.total_pairs}")
print(f"Contraindicated pair   : {report.has_contraindicated}")
print(f"Polypharmacy flagged   : {report.polypharmacy_flagged}")

# 5. View the per-pair breakdown as a DataFrame
df = report_to_dataframe(report)
print(df[["drug_a", "drug_b", "severity", "pair_score"]])
```

### Sample Output

```
Total DDI burden score : 49.0
Risk level             : critical
Interacting pairs      : 5 / 10
Contraindicated pair   : False
Polypharmacy flagged   : True

          drug_a      drug_b severity  pair_score
0     simvastatin  amiodarone    major        14.0
1        warfarin     aspirin    major         7.0
2        warfarin  fluconazole   major         7.0
3        digoxin   amiodarone    major         7.0
...
```

### Custom Thresholds and Weights

```python
# Override severity weights and risk thresholds
custom_weights = {"minor": 1.0, "moderate": 2.0, "major": 5.0, "contraindicated": 20.0}
report = score_regimen(
    patient_meds, db,
    severity_weights=custom_weights,
    high_threshold=15.0,
    critical_threshold=40.0,
    polypharmacy_count=4,
)
```

### Polypharmacy Risk Scorer API

| Function / Class | Description |
|---|---|
| `score_regimen(drugs, db, ...)` | Score a medication list; returns a frozen `RiskReport` |
| `report_to_dataframe(report)` | Flatten a `RiskReport` to a pair-score DataFrame |
| `RiskReport` | Frozen dataclass with score, risk level, and pair details |
| `PairScore` | Frozen dataclass with per-pair severity, weight, and score |
| `SEVERITY_WEIGHTS` | Default weight mapping (`minor=1` → `contraindicated=15`) |
| `HIGH_BURDEN_THRESHOLD` | Score threshold for `"high"` risk (default: 10.0) |
| `CRITICAL_BURDEN_THRESHOLD` | Score threshold for `"critical"` risk (default: 25.0) |
| `POLYPHARMACY_DRUG_COUNT` | Minimum drugs for polypharmacy flag (default: 5) |

---

## New: CYP450 Enzyme-Mediated Interaction Tagger

The `src/cyp450_tagger.py` module scans the `mechanism` text of every DDI
record and labels each row with the cytochrome P450 isoenzymes (and related
transporters such as P-glycoprotein) implicated in the interaction.  Use it
to filter pharmacokinetic interactions by metabolic pathway or to generate
clinician-facing summaries grouped by enzyme.

> **Disclaimer:** illustrative/educational data only — not clinical advice
> and not endorsed by the FDA.

### How It Works

1. Each ``mechanism`` string is matched against a curated catalogue of
   known enzymes (CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, P-glycoprotein,
   and others) using word-boundary-aware, case-insensitive regex.
2. The mechanism is also classified as ``"inhibitor"``, ``"inducer"``, or
   ``"unspecified"`` based on keyword presence.
3. Results are returned as immutable `EnzymeTag` tuples or appended to the
   DataFrame as ``cyp_enzymes`` and ``cyp_role`` columns.

### Step-by-Step Usage

```python
from src import (
    DrugInteractionDB,
    extract_enzymes,
    filter_by_enzyme,
    summarise_by_enzyme,
    tag_interactions,
)

# 1. Load the interaction database
db = DrugInteractionDB()
db.load_data("demo/sample_data.csv")

# 2. Tag every row with CYP enzyme involvement
tagged = tag_interactions(db._df)
print(tagged[["drug_a", "drug_b", "cyp_enzymes", "cyp_role"]].head())

# 3. Aggregate counts by enzyme
summary = summarise_by_enzyme(db._df)
print(summary)

# 4. Find every CYP3A4-inhibitor-mediated interaction
cyp3a4_inhib = filter_by_enzyme(db, "CYP3A4", role="inhibitor")
print(cyp3a4_inhib[["drug_a", "drug_b", "severity", "recommendation"]])

# 5. Tag a single mechanism string ad-hoc
tags = extract_enzymes("Amiodarone inhibits CYP3A4 and CYP2C9")
for t in tags:
    print(f"{t.enzyme}: {t.role}")
```

### Sample Output

```
       drug_a       drug_b   cyp_enzymes   cyp_role
0    warfarin      aspirin                         
1    warfarin  fluconazole        CYP2C9  inhibitor
2 simvastatin   amiodarone CYP3A4, CYP2C9 inhibitor
3    metformin contrast_iodine                     
4 clopidogrel    omeprazole       CYP2C19 inhibitor

  enzyme  inhibitor_count  inducer_count  unspecified_count  total
0 CYP3A4                4              1                  0      5
1 CYP2C9                3              0                  0      3
...
```

### CYP450 Tagger API

| Function / Class | Description |
|---|---|
| `extract_enzymes(mechanism)` | Parse a single mechanism string; returns a tuple of `EnzymeTag` |
| `tag_interactions(df)` | Append `cyp_enzymes` and `cyp_role` columns to a DDI DataFrame |
| `summarise_by_enzyme(df)` | Aggregate inhibitor / inducer / unspecified counts per enzyme |
| `filter_by_enzyme(db, enzyme, role=None)` | Return rows mentioning *enzyme*, optionally restricted by role |
| `EnzymeTag` | Frozen dataclass with `enzyme` and `role` fields |
| `KNOWN_ENZYMES` | Tuple of recognised CYP450 isoforms and transporters |

---

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
