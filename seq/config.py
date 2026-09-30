"""Parts, lengths, paths and the slot scheme the seq/ notebooks share."""

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "out"

PARTS = HERE / "long_parts_fixed.fasta"
BARCODE_MAP = HERE / "barcode_to_circuit_id.csv"

# TODO - To run the genotyping download the data from https://zenodo.org/records/18378852
SUPPLEMENTARY_DIR = HERE.parent / "data/masters_thesis_supplementary"
# Consensus reads, from the Zenodo supplementary data
FASTQ_DIR = SUPPLEMENTARY_DIR / "data/sequencing/2025-01-07_ML_Clone_Validation"

RAW = OUT_DIR / "mseg_raw_michiel.tsv"
FULL_CSV = OUT_DIR / "michiel_lib_genotype_full.csv"
FILTERED_CSV = OUT_DIR / "michiel_lib_genotype_filtered.csv"
REMOVED_CSV = OUT_DIR / "michiel_lib_genotype_removed.csv"
EMPTY_CSV = OUT_DIR / "michiel_lib_genotype_empty.csv"

# Insulators, RBSs and vector padding; dropped before the slot check
STRUCTURAL = {
    "RiboJ", "RiboJ10", "RiboJ51", "RiboJ57",
    "RBS_S1", "RBS_S4", "RBS_P1", "B1_RBS",
    "start_pad", "end_pad",
}

VARIANT = re.compile(r"^(tSrpR|tPhlF)\d+$")


def read_fasta(path):
    """{name: sequence}, keyed on the identifier up to the first space, as mseg keys it."""
    sequences, name, chunks = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


# every record, structural elements and anchors included
FASTA_BP = {name: len(sequence) for name, sequence in read_fasta(PARTS).items()}

# Circuit parts only
PART_BP = {VARIANT.sub(r"\1", name): bp for name, bp in FASTA_BP.items()
           if name not in STRUCTURAL and name not in ("start", "end")}

# Slot options, from michiel_part_library in circmodel/scripts/model.py.
# A circuit is exactly a-b-c-a-b-c; GFP is appended downstream and is not part of the string.
SLOTS = {
    "a": [
        ["AraC", "PhlF"],
        ["pTac", "SrpR"],
    ],
    "b": [
        ["pSrpR", "pTra", "pBAD"],
        ["pSrpR", "pTra"],
        # cryptic site: 
        ["pTac", "tSrpR"],
        ["pTet"],
        ["pTet", "pTac", "tSrpR"],
        ["pTet", "pTac", "tSrpR", "BM3R1"],
        ["Spacer3"],
    ],
    "c": [
        ["BM3R1"],
        ["pTac", "tPhlF"],
        ["pTet", "BM3R1", "pBM3R1"],
        ["pTet", "BM3R1", "pBM3R1", "pBAD"],
        ["pTet"],
        ["pTet", "tPhlF"],
        ["pTra"],
        ["Spacer1"],
        ["Spacer2"],
        ["TraRW"],
    ],
}

# Slot-b option with a cryptic site
CRYPTIC_OPTION = ["pTac", "tSrpR"]
assert CRYPTIC_OPTION in SLOTS["b"]
