"""
Realistic drug-interaction sample data generator.

Run this script to regenerate ``demo/sample_data.csv``:

    python src/data_generator.py

The generator produces pharmacologically accurate drug-drug interaction
records rather than random numeric data, making the output immediately
usable for testing and demonstration.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Any

import pandas as pd


# ---------------------------------------------------------------------------
# Column schema expected by DrugInteractionDB
# ---------------------------------------------------------------------------

COLUMNS: tuple = (
    "drug_a",
    "drug_b",
    "interaction_type",
    "severity",
    "mechanism",
    "clinical_effect",
    "recommendation",
    "evidence_level",
    "source",
)


# ---------------------------------------------------------------------------
# Curated realistic drug-drug interactions
# ---------------------------------------------------------------------------

_INTERACTIONS: List[Dict[str, str]] = [
    {
        "drug_a": "Warfarin",
        "drug_b": "Aspirin",
        "interaction_type": "pharmacodynamic",
        "severity": "major",
        "mechanism": "Additive anticoagulant and antiplatelet effects",
        "clinical_effect": "Increased risk of bleeding including gastrointestinal hemorrhage",
        "recommendation": "Avoid combination; if unavoidable monitor INR closely",
        "evidence_level": "A",
        "source": "Hansten & Horn Drug Interactions Analysis",
    },
    {
        "drug_a": "Warfarin",
        "drug_b": "Fluconazole",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Fluconazole inhibits CYP2C9 reducing warfarin metabolism",
        "clinical_effect": "Elevated warfarin plasma levels leading to INR increase and bleeding risk",
        "recommendation": "Monitor INR closely; reduce warfarin dose by 25-50% as needed",
        "evidence_level": "A",
        "source": "Clinical Pharmacokinetics 2001",
    },
    {
        "drug_a": "Simvastatin",
        "drug_b": "Amiodarone",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Amiodarone inhibits CYP3A4 and CYP2C9 reducing statin clearance",
        "clinical_effect": "Markedly elevated statin levels increasing risk of myopathy and rhabdomyolysis",
        "recommendation": "Do not exceed simvastatin 20 mg/day; consider alternative statin",
        "evidence_level": "A",
        "source": "FDA Drug Safety Communication 2011",
    },
    {
        "drug_a": "Metformin",
        "drug_b": "Contrast_Iodine",
        "interaction_type": "pharmacodynamic",
        "severity": "major",
        "mechanism": "Contrast-induced nephropathy reduces metformin excretion",
        "clinical_effect": "Accumulation of metformin leading to lactic acidosis risk",
        "recommendation": "Hold metformin 48 hours before and after iodinated contrast administration",
        "evidence_level": "B",
        "source": "ACR Manual on Contrast Media 2022",
    },
    {
        "drug_a": "Clopidogrel",
        "drug_b": "Omeprazole",
        "interaction_type": "pharmacokinetic",
        "severity": "moderate",
        "mechanism": "Omeprazole inhibits CYP2C19 reducing clopidogrel activation",
        "clinical_effect": "Reduced antiplatelet effect increasing cardiovascular event risk",
        "recommendation": "Use pantoprazole instead; avoid omeprazole with clopidogrel",
        "evidence_level": "A",
        "source": "NEJM 2010 COGENT Trial",
    },
    {
        "drug_a": "Ciprofloxacin",
        "drug_b": "Theophylline",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Ciprofloxacin inhibits CYP1A2 reducing theophylline clearance",
        "clinical_effect": "Theophylline toxicity including nausea, vomiting, seizures, and arrhythmias",
        "recommendation": "Reduce theophylline dose 30-50%; monitor serum levels closely",
        "evidence_level": "A",
        "source": "Clinical Pharmacokinetics 1992",
    },
    {
        "drug_a": "Lithium",
        "drug_b": "Ibuprofen",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "NSAIDs reduce renal prostaglandin synthesis decreasing lithium clearance",
        "clinical_effect": "Elevated lithium plasma concentrations leading to toxicity",
        "recommendation": "Avoid NSAIDs in lithium patients; use acetaminophen instead; monitor levels",
        "evidence_level": "A",
        "source": "Pharmacotherapy 2003",
    },
    {
        "drug_a": "Digoxin",
        "drug_b": "Amiodarone",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Amiodarone inhibits P-glycoprotein and renal tubular secretion of digoxin",
        "clinical_effect": "Digoxin toxicity with bradycardia, nausea, visual disturbances, and arrhythmias",
        "recommendation": "Reduce digoxin dose by 50% at initiation; monitor serum levels and ECG",
        "evidence_level": "A",
        "source": "Journal of Clinical Pharmacology 1985",
    },
    {
        "drug_a": "SSRIs",
        "drug_b": "MAOIs",
        "interaction_type": "pharmacodynamic",
        "severity": "contraindicated",
        "mechanism": "Additive serotonergic activity causing serotonin syndrome",
        "clinical_effect": "Life-threatening serotonin syndrome with hyperthermia, rigidity, and seizures",
        "recommendation": "Contraindicated; allow 14-day washout between agents",
        "evidence_level": "A",
        "source": "FDA Label and Case Reports",
    },
    {
        "drug_a": "Methotrexate",
        "drug_b": "Trimethoprim",
        "interaction_type": "pharmacodynamic",
        "severity": "major",
        "mechanism": "Additive folate antagonism increasing toxicity risk",
        "clinical_effect": "Severe bone marrow suppression, megaloblastic anemia, and mucositis",
        "recommendation": "Avoid combination; if needed monitor CBC and folate levels closely",
        "evidence_level": "B",
        "source": "British Journal of Clinical Pharmacology 2007",
    },
    {
        "drug_a": "Sildenafil",
        "drug_b": "Nitrates",
        "interaction_type": "pharmacodynamic",
        "severity": "contraindicated",
        "mechanism": "Additive cGMP-mediated vasodilation causing profound hypotension",
        "clinical_effect": "Severe and potentially fatal hypotension requiring emergency intervention",
        "recommendation": "Absolute contraindication; do not administer within 24 hours of nitrate use",
        "evidence_level": "A",
        "source": "FDA Contraindication Label",
    },
    {
        "drug_a": "Carbamazepine",
        "drug_b": "Valproate",
        "interaction_type": "pharmacokinetic",
        "severity": "moderate",
        "mechanism": "Mutual enzyme induction and inhibition altering both drug levels",
        "clinical_effect": "Carbamazepine toxicity and reduced valproate efficacy with unpredictable levels",
        "recommendation": "Monitor serum levels of both drugs and adjust doses accordingly",
        "evidence_level": "A",
        "source": "Epilepsia 2006",
    },
    {
        "drug_a": "Rifampicin",
        "drug_b": "Oral_Contraceptives",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Rifampicin induces CYP3A4 increasing contraceptive hormone metabolism",
        "clinical_effect": "Reduced contraceptive efficacy leading to unintended pregnancy",
        "recommendation": "Use alternative or additional non-hormonal contraception during rifampicin therapy",
        "evidence_level": "A",
        "source": "Obstetrics and Gynecology 2011",
    },
    {
        "drug_a": "ACE_Inhibitors",
        "drug_b": "Potassium_Sparing_Diuretics",
        "interaction_type": "pharmacodynamic",
        "severity": "moderate",
        "mechanism": "Additive effect on renal potassium retention",
        "clinical_effect": "Hyperkalemia with risk of cardiac arrhythmias especially in renal impairment",
        "recommendation": "Monitor serum potassium levels; avoid in patients with renal impairment",
        "evidence_level": "B",
        "source": "Journal of Hypertension 2004",
    },
    {
        "drug_a": "Fluoxetine",
        "drug_b": "Tamoxifen",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Fluoxetine strongly inhibits CYP2D6 impairing tamoxifen activation to endoxifen",
        "clinical_effect": "Reduced tamoxifen efficacy potentially increasing breast cancer recurrence risk",
        "recommendation": "Avoid combination; use venlafaxine or gabapentin for hot flashes instead",
        "evidence_level": "A",
        "source": "JNCI 2010",
    },
    {
        "drug_a": "Linezolid",
        "drug_b": "Pseudoephedrine",
        "interaction_type": "pharmacodynamic",
        "severity": "major",
        "mechanism": "Additive adrenergic effects due to MAO inhibition by linezolid",
        "clinical_effect": "Hypertensive crisis with severe headache and risk of stroke",
        "recommendation": "Avoid sympathomimetics during linezolid therapy and for 2 weeks after",
        "evidence_level": "B",
        "source": "Antimicrobial Agents and Chemotherapy 2004",
    },
    {
        "drug_a": "Tacrolimus",
        "drug_b": "Fluconazole",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Fluconazole inhibits CYP3A4 and P-glycoprotein reducing tacrolimus metabolism",
        "clinical_effect": "Tacrolimus toxicity with nephrotoxicity, neurotoxicity, and immunosuppression excess",
        "recommendation": "Reduce tacrolimus dose by 50%; monitor trough levels every 2-3 days",
        "evidence_level": "A",
        "source": "Transplantation 2002",
    },
    {
        "drug_a": "Methadone",
        "drug_b": "Fluconazole",
        "interaction_type": "pharmacokinetic",
        "severity": "moderate",
        "mechanism": "Fluconazole inhibits CYP3A4 and CYP2C19 reducing methadone clearance",
        "clinical_effect": "Elevated methadone levels with sedation, respiratory depression, and QTc prolongation",
        "recommendation": "Monitor for methadone toxicity; reduce dose if QTc prolongation observed",
        "evidence_level": "B",
        "source": "Journal of Pain and Symptom Management 2008",
    },
    {
        "drug_a": "Allopurinol",
        "drug_b": "Azathioprine",
        "interaction_type": "pharmacokinetic",
        "severity": "major",
        "mechanism": "Allopurinol inhibits xanthine oxidase essential for azathioprine metabolism",
        "clinical_effect": "Severe accumulation of azathioprine metabolites causing bone marrow suppression",
        "recommendation": "Reduce azathioprine dose by 75%; monitor CBC weekly for first month",
        "evidence_level": "A",
        "source": "Annals of Rheumatic Diseases 1966",
    },
    {
        "drug_a": "Phenytoin",
        "drug_b": "Folic_Acid",
        "interaction_type": "pharmacokinetic",
        "severity": "moderate",
        "mechanism": "Folic acid may increase phenytoin metabolism and reduce absorption",
        "clinical_effect": "Lowered phenytoin serum levels risking breakthrough seizures",
        "recommendation": "Monitor phenytoin levels when starting or stopping folate supplementation",
        "evidence_level": "B",
        "source": "Epilepsia 1999",
    },
]


def generate_sample(
    n: int = 20,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Return a DataFrame of curated, realistic drug-drug interaction records.

    Unlike a purely random generator, this function returns a fixed set of
    pharmacologically accurate interactions drawn from published literature.
    The *n* and *seed* parameters are accepted for API compatibility but only
    affect how many rows are returned when *n* is smaller than the full set.

    Args:
        n: Maximum number of rows to return.  When *n* is greater than or
           equal to the number of available interactions all rows are
           returned.  Defaults to 20.
        seed: Unused — retained for API compatibility with callers that
              pass a random seed.  The curated dataset is deterministic.

    Returns:
        DataFrame with columns matching :data:`COLUMNS`.
    """
    rows = _INTERACTIONS[:n]
    return pd.DataFrame(rows, columns=list(COLUMNS))


if __name__ == "__main__":
    Path("demo").mkdir(exist_ok=True)
    df = generate_sample(20)
    out_path = "demo/sample_data.csv"
    df.to_csv(out_path, index=False)
    print(f"Generated {len(df)} records -> {out_path}")
    print(df[["drug_a", "drug_b", "severity"]].to_string(index=False))
    print(f"\nShape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
