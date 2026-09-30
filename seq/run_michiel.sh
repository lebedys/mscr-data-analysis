#!/usr/bin/env bash
# Genotype the Michiel clone-validation consensus reads.
#
# One mseg call per barcode: each .final.fastq holds a single assembled consensus, not a read set.
# -d requires both anchors, which is what drops the reads that never assembled across the junction;
# --circular is inert without it, because the rotation branch sits inside the -d arm.
#
# MSEG must be the patched binary: it emits a third `coverage` column, and the notebook keeps only
# raw lines with exactly three fields. An unpatched build emits two columns and so yields an empty
# report rather than an error.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MSEG="${MSEG:-$HOME/mscr/conrad/mseg/target/release/mseg}"
FASTQ_DIR="${FASTQ_DIR:-$HERE/../data/masters_thesis_supplementary/data/sequencing/2025-01-07_ML_Clone_Validation}"
PARTS="${PARTS:-$HERE/long_parts_fixed.fasta}"
OUT="${OUT:-$HERE/out/mseg_raw_michiel.tsv}"
LOG="${LOG:-$HERE/out/mseg_michiel.log}"

mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")"
: > "$OUT"
: > "$LOG"
TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT

for fastq in "$FASTQ_DIR"/barcode*.final.fastq; do
    "$MSEG" \
        --segments "$PARTS" \
        --sequences "$fastq" \
        --min-seq-len 200 \
        --max-seq-len 50000 \
        --min-norm-score 1.3 \
        --circular \
        -d \
        --output "$TMP" 2>> "$LOG"
    cat "$TMP" >> "$OUT"
done

echo "$OUT: $(wc -l < "$OUT") genotypes from $(ls "$FASTQ_DIR"/barcode*.final.fastq | wc -l) barcodes"
