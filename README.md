# KRAS Variant Classifier

A small Python script that classifies single-nucleotide variants in **KRAS** by their effect on the protein: **Silent**, **Missense**, **Nonsense**, **Stop_loss** or **Start_loss**.

It downloads the coding sequence (CDS) of the RefSeq transcript `NM_004985.5` from NCBI, applies each variant to the affected codon, translates the original and altered codons, and compares the resulting amino acids.

## How it works

1. Fetches the GenBank record for `NM_004985.5` through NCBI Entrez and extracts the CDS.
2. Reads the variant table and parses the coding position from each `variant_name` (e.g. `c.35G>A` gives position 35).
3. **Strand validation:** KRAS is on the minus strand, so ClinVar-style `ref`/`alt` bases are often the complement of the transcript bases. Each variant's `ref` is compared with the CDS base and labelled:
   - `DIRECT_MATCH`: `ref` equals the CDS base.
   - `COMPLEMENT_MATCH`: the complement of `ref` equals the CDS base, so `alt` is complemented before use.
   - `MISMATCH`: neither matches, so the variant is not classified.
4. Substitutes the base into the affected codon, translates both codons and classifies the change.

Variants that cannot be resolved are kept in the output as `Unresolved`, with the reason in `match_status` (`NO_POSITION_PARSED`, `OUT_OF_RANGE`, `NOT_SNV` or `MISMATCH`).

## Requirements

- Python 3.8+
- [Biopython](https://biopython.org/) and [pandas](https://pandas.pydata.org/)

```bash
pip install -r requirements.txt
```

## Usage

Before running, open `classify_kras_variants.py` and edit two lines near the top and in the NCBI section:

1. **Input file path:** set `VARIANTS_CSV` to the path of your variants CSV. To run the bundled demo, use:
   ```python
   VARIANTS_CSV = "kras_variants.csv"
   ```
   For a file elsewhere on your machine, use a raw string so Windows backslashes work:
   ```python
   VARIANTS_CSV = r"C:\path\to\your_variants.csv"
   ```
2. **Entrez email:** NCBI asks every Entrez user to identify themselves, so replace the placeholder with your own email address:
   ```python
   Entrez.email = "you@example.com"
   ```

Then run the script from the repo folder:

```bash
python classify_kras_variants.py
```

The results are printed to the terminal and saved to `kras_variants_classified.csv` in the current folder. To change the output name, edit `OUTPUT_CSV`.

> **Note:** if you fork or contribute to this repo, don't commit your real email address. Keep the placeholder in the committed version.

## Input format

A CSV with at least these columns:

| Column         | Description                                                                |
|----------------|----------------------------------------------------------------------------|
| `variant_name` | HGVS-style name containing a coding change, e.g. `NM_004985.5(KRAS):c.35G>A` |
| `ref`          | Reference base (single nucleotide)                                         |
| `alt`          | Alternate base (single nucleotide)                                         |

Any other columns (such as `position`, `gene`, `clinical_significance`) are carried through to the output unchanged. The `position` column is not used by the script; the coding position is parsed from `variant_name`.

A demo file, `kras_variants.csv`, is included. It contains KRAS variants from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) with GRCh38 positions.

## Output

The input columns plus:

| Column            | Description                                                              |
|-------------------|--------------------------------------------------------------------------|
| `c_positions`     | Coding position parsed from `variant_name`                               |
| `match_status`    | Strand validation result (see above)                                     |
| `original_aa`     | Reference amino acid (single-letter code, `*` = stop)                    |
| `altered_aa`      | Amino acid after the variant                                             |
| `classification`  | `Silent`, `Missense`, `Nonsense`, `Stop_loss`, `Start_loss` or `Unresolved` |

Expected results for a few of the demo variants:

| Variant | Protein change | Classification |
|---------|----------------|----------------|
| `c.34G>T` | Gly12Cys (G12C) | Missense |
| `c.35G>A` | Gly12Asp (G12D) | Missense |
| `c.35G>T` | Gly12Val (G12V) | Missense |
| `c.183A>C` | Gln61His (Q61H) | Missense |

`c.35G>A` giving G12D and `c.34G>T` giving G12C are good quick checks that the script works.

## Limitations

- Handles single-nucleotide substitutions only. Insertions, deletions and multi-base changes are reported as `NOT_SNV`.
- Handles coding (`c.`) positions only. Intronic and UTR notations such as `c.123+2T>C` or `c.-5A>G` are reported as `NO_POSITION_PARSED`.
- Classification is based on the translated codon only; it does not predict splicing effects or clinical impact.
- Requires an internet connection to fetch the transcript from NCBI.

## Data source

Demo variants come from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (NCBI), which is in the public domain. Transcript sequence is retrieved from NCBI RefSeq via Entrez.

## Author

Vanssha Kaitha
