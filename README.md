# Combinatorial Synthesis of Diverse Genetic Circuit Libraries (Data Analysis)

Notebooks under the `nb/` directory contain the characterisation pipeline described in Chapter 3 and applied in Chapter 4.
Global configuration parameters for these noteboks are defined in `nb/analyse_pr_config.py`.

Genotyping scripts and notebooks are included in `seq/`.
Global configuration parameters for these scripts/notebooks are defined in `seq/config.py`.

## Plate Reader Analysis (`nb/`)

### `analyse_pr_biomass.ipynb`
Growth rate of every well by three methods (polyfit, Savitzky-Golay, moving window), and the timepoint each samples at.
- input: plate-reader files for CD1, CD2, AL1 (see [Data](#data))
- output: `out/analyse_pr_biomass.csv` (growth rate, sampled timepoint and biomass per group, replicate, condition, and method)
- output: `out/analyse_pr_biomass_{polyfit,sg,mwindow}.pdf` (growth rate and ln(B) per group, one PDF per method)
- output: `out/analyse_pr_biomass_bg_correction.pdf` (raw against blank-corrected ln(OD) per group)

### `analyse_pr_biomass_naive.ipynb`
Compares the SG sampled timepoint against naive sampling by time or OD thresholds.
- input: AL1 plate-reader file (see [Data](#data))
- input: `out/analyse_pr_biomass.csv` (from `analyse_pr_biomass.ipynb`)
- output: `out/analyse_pr_biomass_naive.pdf` (figure of failure modes of naive methods)

### `analyse_pr_fl.ipynb`
Fluorescence over time for every group.
- input: plate-reader files for CD1, CD2, AL1
- input: `out/analyse_pr_biomass.csv` (from `analyse_pr_biomass.ipynb`)
- output: `out/analyse_pr_fl.pdf` (F/OD, F and ln(OD) per group)

### `analyse_pr_autof_pL2.ipynb`
Builds pL2 autofluorescence calibration curve.
- input: plate-reader files for CD1, CD2, AL1
- input: `out/analyse_pr_biomass.csv` (from `analyse_pr_biomass.ipynb`)
- output: n/a

### `analyse_pr_synthrate.ipynb`
Synthesis rate and RPU of every group at its sampled timepoint.
- input: plate-reader files for CD1, CD2, AL1
- input: `out/analyse_pr_biomass.csv` (from `analyse_pr_biomass.ipynb`)
- output: `out/analyse_pr_synthrate.csv` (synthesis rate and RPU per group, replicate, condition, and method)
- output: `out/analyse_pr_synthrate_plots.pdf` (synthrate, dG/dt, F and ln(B) per group and replicate)
- output: `out/analyse_pr_synthrate_bars.pdf` (RPU per condition, six groups per page)
- output: `out/analyse_pr_synthrate_bars/*.png` (the same chart, one PNG per group)

### `analyse_pr_synthrate_scatter.ipynb`
How much the growth-rate method changes the synthesis rate and RPU.
- input: `out/analyse_pr_synthrate.csv` (from `analyse_pr_synthrate.ipynb`)
- output: n/a

### `analyse_pr_classify.ipynb`
Classifies circuit logic gate by KL divergence and by Cello score.
- input: `out/analyse_pr_synthrate.csv` (from `analyse_pr_synthrate.ipynb`)
- input: `seq/flow_invalid.txt` (circuits with multimodal flow distributions)
- output: `out/analyse_pr_classify.csv` (gate calls, per-condition RPU and scores per group)
- output: `out/analyse_pr_classify.pdf` (RPU bars and both rankings, one page per group)
- output: `out/analyse_pr_classify/*.png` (the same chart, one PNG per group)

### `analyse_flow.ipynb`
MEFL distributions of every circuit from flow cytometry.
- input: CD1, CD2 flow cytometry FCS files and `sample_to_circuit.csv` (see [Data](#data))
- output: `out/analyse_flow/*.png` (MEFL density per replicate, one PNG per group)

### `titration_75/analyse_pr_titration.ipynb`
Growth, RPU and dose response of Circuit 75 across an aTc × IPTG grid.
- input: `data/plate_reader/20251104_titration_iptg_atc_p75-CD-ML3/atc_iptg_titration_20251104.txt`
- output: `titration_75/out/analyse_pr_titration.csv` (RPU per aTc and IPTG level)
- output: `titration_75/out/titration_{blanks,controls,growth,calibration,bars,dose_response}.pdf` (one figure per pipeline step)

## Genotyping (`seq/`)

### `run_michiel.sh`
Calls parts in each Nanopore consensus read with mseg.
- input: consensus FASTQ reads (see [Data](#data))
- input: `long_parts_fixed.fasta` (the parts vocabulary)
- output: `out/mseg_raw_michiel.tsv` (called parts per barcode)

### `michiel_lib_genotype.ipynb`
Builds the genotype table and filters it on read coverage.
- input: `out/mseg_raw_michiel.tsv` (from `run_michiel.sh`)
- input: consensus FASTQ reads, `barcode_to_circuit_id.csv`
- input: `nb/out/analyse_pr_classify.csv` (from `analyse_pr_classify.ipynb`)
- output: `out/michiel_lib_genotype_full.csv` (every barcode)
- output: `out/michiel_lib_genotype_filtered.csv` (the barcodes that pass every filter)
- output: `out/michiel_lib_genotype_removed.csv` (genotyped barcodes a filter removed)
- output: `out/michiel_lib_genotype_empty.csv` (barcodes with no genotype)

### `analyse_parts.ipynb`
Part and slot usage across the sequenced circuits.
- input: `out/michiel_lib_genotype_filtered.csv`, `out/michiel_lib_genotype_removed.csv` (from `michiel_lib_genotype.ipynb`)
- output: `out/part_histo.png`, `out/part_scatter.png` (part occurrences, and against length)
- output: `out/composites_histo.png`, `out/composites_scatter.png` (slot option usage, and against length)
- output: `out/part_availability.png`, `out/part_availability_cryptic.png` (use against availability, without and with the cryptic option)
- output: `out/part_histo_removed.png`, `out/part_histo_invalid.png` (part occurrences in the excluded circuits)

### `plot_sbol.ipynb`
SBOL diagram of every filtered circuit.
- input: `out/michiel_lib_genotype_filtered.csv` (from `michiel_lib_genotype.ipynb`)
- input: `nb/out/analyse_pr_classify/*.png`, `nb/out/analyse_pr_synthrate_bars/*.png`, `nb/out/analyse_flow/*.png`
- output: `out/sbol/*.png` (one diagram per circuit)
- output: `out/classify_sbol/*.png`, `out/bars_sbol/*.png`, `out/classify_flow_sbol/*.png` (diagram under each chart)
- output: `out/classify_sbol.pdf`, `out/classify_flow_sbol.pdf` (the classify composites, one page per circuit)

## Data

- AL1 plate reader: `data/plate_reader/2025_10_21_ML3_384/`
- Titration plate reader: `data/plate_reader/20251104_titration_iptg_atc_p75-CD-ML3/`
- Titration plate layout: `data/protocols/atc_iptg_titration.xlsx`
- CD1, CD2 plate reader, CD1, CD2 flow cytometry and the sequencing reads:
  [Zenodo record 18378852](https://zenodo.org/records/18378852). Unzip `masters_thesis_supplementary.zip`
  into `data/`, keeping its layout.
